/**
 * DICOM File Parser & 2D Image Loader
 * HydroMorph — hydrocephalus morphometrics pipeline
 * Supports: .dcm (DICOM), .png, .jpg, .jpeg, .bmp (single-slice images)
 *
 * Requires: dicomParser (https://unpkg.com/dicom-parser)
 */

const DicomReader = (() => {

  /**
   * Parse a single DICOM file (.dcm) from ArrayBuffer.
   * Returns { pixelData: TypedArray, rows, cols, sliceLocation, instanceNumber,
   *           rescaleSlope, rescaleIntercept, spacingX, spacingY, sliceThickness }
   */
  function parseFile(arrayBuffer) {
    const byteArray = new Uint8Array(arrayBuffer);
    const dataSet = dicomParser.parseDicom(byteArray);

    // Extract metadata
    const rows    = dataSet.uint16('x00280010') || 0;
    const cols    = dataSet.uint16('x00280011') || 0;
    const bitsAllocated       = dataSet.uint16('x00280100') || 16;
    const pixelRepresentation = dataSet.uint16('x00280103') || 0; // 0 = unsigned, 1 = signed
    const rescaleSlope        = parseFloat(dataSet.string('x00281053') || '1');
    const rescaleIntercept    = parseFloat(dataSet.string('x00281052') || '0');
    const sliceLocation       = parseFloat(dataSet.string('x00201041') || '0');
    const instanceNumber      = parseInt(dataSet.string('x00200013') || '0', 10);
    const pixelSpacingStr     = dataSet.string('x00280030') || '1\\1';
    const sliceThickness      = parseFloat(dataSet.string('x00180050') || '1');

    if (rows === 0 || cols === 0) {
      throw new Error('DICOM file has invalid dimensions (rows/cols = 0).');
    }

    // Parse pixel spacing (row spacing \ col spacing → spacingY, spacingX in mm)
    const parts = pixelSpacingStr.split('\\').map(Number);
    const spacingY = (parts[0] > 0 && isFinite(parts[0])) ? parts[0] : 1;
    const spacingX = (parts[1] > 0 && isFinite(parts[1])) ? parts[1] : 1;

    // Get pixel data element
    const pixelDataElement = dataSet.elements.x7fe00010;
    if (!pixelDataElement) {
      throw new Error('No pixel data (7FE0,0010) found in DICOM file.');
    }

    const numPixels = rows * cols;
    const pixelData = extractPixelData(
      byteArray,
      pixelDataElement,
      numPixels,
      bitsAllocated,
      pixelRepresentation
    );

    return {
      pixelData,
      rows,
      cols,
      sliceLocation,
      instanceNumber,
      rescaleSlope:    isFinite(rescaleSlope) ? rescaleSlope : 1,
      rescaleIntercept: isFinite(rescaleIntercept) ? rescaleIntercept : 0,
      spacingX,
      spacingY,
      sliceThickness: (sliceThickness > 0 && isFinite(sliceThickness)) ? sliceThickness : 1,
    };
  }

  /**
   * Extract pixel data from a DICOM dataset, handling 8-bit, 16-bit unsigned,
   * and 16-bit signed representations.
   */
  function extractPixelData(byteArray, element, numPixels, bitsAllocated, pixelRepresentation) {
    const offset = element.dataOffset;
    const length = element.length;

    if (bitsAllocated === 8) {
      // 8-bit — always unsigned in practice
      const result = new Float32Array(numPixels);
      for (let i = 0; i < numPixels && i < length; i++) {
        result[i] = byteArray[offset + i];
      }
      return result;
    }

    if (bitsAllocated === 16) {
      const result = new Float32Array(numPixels);
      const maxI = Math.min(numPixels, Math.floor(length / 2));
      if (pixelRepresentation === 1) {
        // Signed 16-bit (Int16)
        for (let i = 0; i < maxI; i++) {
          const lo = byteArray[offset + i * 2];
          const hi = byteArray[offset + i * 2 + 1];
          let val = (hi << 8) | lo; // little-endian
          if (val >= 0x8000) val -= 0x10000; // sign extend
          result[i] = val;
        }
      } else {
        // Unsigned 16-bit (Uint16)
        for (let i = 0; i < maxI; i++) {
          const lo = byteArray[offset + i * 2];
          const hi = byteArray[offset + i * 2 + 1];
          result[i] = (hi << 8) | lo; // little-endian
        }
      }
      return result;
    }

    if (bitsAllocated === 32) {
      // 32-bit — treat as float
      const buf = byteArray.buffer.slice(
        byteArray.byteOffset + offset,
        byteArray.byteOffset + offset + numPixels * 4
      );
      return new Float32Array(buf);
    }

    throw new Error(`Unsupported BitsAllocated: ${bitsAllocated}`);
  }

  /**
   * Parse multiple DICOM files (a series) into a 3D volume.
   * Sorts slices by SliceLocation, then InstanceNumber as fallback.
   * Returns { shape, spacing, affine, data, header } — same format as NiftiReader.parse()
   */
  async function parseSeries(arrayBuffers) {
    if (!arrayBuffers || arrayBuffers.length === 0) {
      throw new Error('No DICOM files provided.');
    }

    // Parse each file
    const slices = [];
    for (let i = 0; i < arrayBuffers.length; i++) {
      try {
        slices.push(parseFile(arrayBuffers[i]));
      } catch (e) {
        console.warn(`Skipping DICOM file ${i}: ${e.message}`);
      }
    }

    if (slices.length === 0) {
      throw new Error('No valid DICOM slices could be parsed.');
    }

    // Sort by sliceLocation, fall back to instanceNumber
    slices.sort((a, b) => {
      const dLoc = a.sliceLocation - b.sliceLocation;
      if (Math.abs(dLoc) > 0.001) return dLoc;
      return a.instanceNumber - b.instanceNumber;
    });

    const rows     = slices[0].rows;
    const cols     = slices[0].cols;
    const numSlices = slices.length;
    const spacingX = slices[0].spacingX || 1;
    const spacingY = slices[0].spacingY || 1;

    // Calculate slice spacing from sorted locations
    let spacingZ = slices[0].sliceThickness || 1;
    if (slices.length > 1) {
      // Use median of consecutive differences for robustness
      const diffs = [];
      for (let i = 1; i < slices.length; i++) {
        const dz = Math.abs(slices[i].sliceLocation - slices[i - 1].sliceLocation);
        if (dz > 0.001) diffs.push(dz);
      }
      if (diffs.length > 0) {
        diffs.sort((a, b) => a - b);
        spacingZ = diffs[Math.floor(diffs.length / 2)]; // median
      }
    }

    // Stack into 3D volume using the same layout as NIfTI reader:
    // index = x + y * cols + z * cols * rows
    const total = cols * rows * numSlices;
    const data  = new Float32Array(total);

    for (let z = 0; z < numSlices; z++) {
      const slice = slices[z];
      const { rescaleSlope: slope, rescaleIntercept: intercept, pixelData } = slice;

      // Make sure this slice has compatible dimensions
      if (slice.rows !== rows || slice.cols !== cols) {
        console.warn(`Slice ${z} has different dimensions (${slice.cols}×${slice.rows} vs ${cols}×${rows}), skipping`);
        continue;
      }

      for (let y = 0; y < rows; y++) {
        for (let x = 0; x < cols; x++) {
          const rawVal = pixelData[y * cols + x];
          // Apply rescale: HU = pixel * slope + intercept
          data[x + y * cols + z * cols * rows] = rawVal * slope + intercept;
        }
      }
    }

    return {
      shape:   [cols, rows, numSlices],
      spacing: [spacingX, spacingY, spacingZ],
      affine: [
        [spacingX, 0, 0, 0],
        [0, spacingY, 0, 0],
        [0, 0, spacingZ, 0],
        [0, 0, 0, 1]
      ],
      data,
      header: {
        ndim:     3,
        datatype: 16,
        bitpix:   32,
        voxOffset: 0,
        sformCode: 0,
        dims:    [cols, rows, numSlices],
        pixdim:  [spacingX, spacingY, spacingZ],
      },
    };
  }

  /**
   * Parse a single 2D image (PNG/JPG/BMP) as a grayscale "single-slice" pseudo-volume.
   * Returns { shape, spacing, affine, data, header, isSingleSlice: true }
   */
  function parseImage(file) {
    return new Promise((resolve, reject) => {
      const objectUrl = URL.createObjectURL(file);
      const img = new Image();

      img.onload = () => {
        try {
          URL.revokeObjectURL(objectUrl);

          const W = img.naturalWidth  || img.width;
          const H = img.naturalHeight || img.height;

          const canvas = document.createElement('canvas');
          canvas.width  = W;
          canvas.height = H;
          const ctx = canvas.getContext('2d');
          ctx.drawImage(img, 0, 0);

          const imageData = ctx.getImageData(0, 0, W, H);
          const total = W * H;
          const data  = new Float32Array(total);

          // Luma-weighted grayscale (0–255 range, acts as pseudo-HU)
          for (let i = 0; i < total; i++) {
            const r = imageData.data[i * 4];
            const g = imageData.data[i * 4 + 1];
            const b = imageData.data[i * 4 + 2];
            data[i] = 0.299 * r + 0.587 * g + 0.114 * b;
          }

          resolve({
            shape:   [W, H, 1],
            spacing: [1, 1, 1],
            affine: [
              [1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1],
            ],
            data,
            header: {
              ndim:     2,
              datatype: 16,
              bitpix:   32,
              voxOffset: 0,
              sformCode: 0,
              dims:    [W, H, 1],
              pixdim:  [1, 1, 1],
            },
            isSingleSlice: true,
          });
        } catch (e) {
          reject(new Error(`Failed to process image pixels: ${e.message}`));
        }
      };

      img.onerror = () => {
        URL.revokeObjectURL(objectUrl);
        reject(new Error('Failed to load image file.'));
      };

      img.src = objectUrl;
    });
  }

  return { parseFile, parseSeries, parseImage };

})();
