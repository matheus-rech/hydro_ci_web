/**
 * MedSAM2 API Client — Gradio HTTP API
 * HydroMorph — hydrocephalus morphometrics pipeline
 *
 * Communicates with a MedSAM2 Gradio server on Hugging Face Spaces (ZeroGPU).
 * Uses the plain-fetch Gradio API pattern:
 *   1. Upload file → POST /gradio_api/upload
 *   2. Submit job  → POST /call/<api_name>
 *   3. Read result → GET  /call/<api_name>/<event_id> (SSE stream)
 *
 * All segmentation is optional — the app falls back to the threshold
 * pipeline when the server is unavailable.
 *
 * Requires: pako (for gzip compression/decompression)
 */

const MedSAMClient = (() => {

  // Default: Hugging Face Space URL
  let serverUrl = 'https://mmrech-medsam2-server.hf.space';

  /** Update the server base URL (strips trailing slash). */
  function setServerUrl(url) {
    serverUrl = (url || '').trim().replace(/\/$/, '');
  }

  /** Return the current server URL. */
  function getServerUrl() {
    return serverUrl;
  }

  // ─── Gradio SSE Reader ────────────────────────────────────────────────────

  /**
   * Read a Gradio SSE stream and return the final data array.
   * Handles events: complete, error, heartbeat, generating.
   */
  async function readGradioSSE(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`SSE fetch failed: ${resp.status}`);
    const reader = resp.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    let eventType = null;

    while (true) {
      const { done, value } = reader.read ? await reader.read() : { done: true };
      if (done) break;
      buffer += decoder.decode(value, { stream: true });

      const lines = buffer.split('\n');
      buffer = lines.pop();

      for (const line of lines) {
        if (line.startsWith('event: ')) {
          eventType = line.slice(7).trim();
        } else if (line.startsWith('data: ') && eventType) {
          const data = JSON.parse(line.slice(6));
          if (eventType === 'complete') {
            reader.cancel();
            return data;
          }
          if (eventType === 'error') {
            reader.cancel();
            throw new Error('Gradio API error: ' + JSON.stringify(data));
          }
          // 'generating' / 'heartbeat' events — ignore and continue
        }
      }
    }
    throw new Error('SSE stream ended without a complete event');
  }

  // ─── File Upload ──────────────────────────────────────────────────────────

  /**
   * Upload a File/Blob to the Gradio server.
   * Returns a Gradio FileData reference object.
   */
  async function uploadFile(blob, filename) {
    const fd = new FormData();
    fd.append('files', blob, filename);

    const resp = await fetch(`${serverUrl}/gradio_api/upload`, {
      method: 'POST',
      body: fd,
    });
    if (!resp.ok) throw new Error(`Upload failed: ${resp.status}`);
    const paths = await resp.json();

    return {
      path: paths[0],
      orig_name: filename,
      meta: { '_type': 'gradio.FileData' },
    };
  }

  // ─── Submit + Read ────────────────────────────────────────────────────────

  /**
   * Call a Gradio API endpoint (POST → SSE pattern).
   * Returns the parsed data array from the 'complete' event.
   */
  async function callGradioApi(apiName, dataArray) {
    const submitResp = await fetch(`${serverUrl}/call/${apiName}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ data: dataArray }),
    });
    if (!submitResp.ok) {
      const txt = await submitResp.text().catch(() => '');
      throw new Error(`Gradio submit failed (${submitResp.status}): ${txt}`);
    }
    const { event_id } = await submitResp.json();
    return readGradioSSE(`${serverUrl}/call/${apiName}/${event_id}`);
  }

  // ─── NPY Serialisation (browser) ─────────────────────────────────────────

  /**
   * Serialise a typed array into a minimal .npy (NumPy v1.0) file in the
   * browser. Supports float32 and uint8.
   *
   * @param {Float32Array|Uint8Array} data - flat array
   * @param {number[]} shape - e.g. [D, H, W]
   * @returns {Blob} .npy file as a Blob
   */
  function createNpyBlob(data, shape) {
    let dtype;
    if (data instanceof Float32Array) dtype = '<f4';
    else if (data instanceof Uint8Array) dtype = '|u1';
    else if (data instanceof Int16Array) dtype = '<i2';
    else throw new Error('Unsupported typed array for npy serialisation');

    const shapeStr = shape.length === 1
      ? `(${shape[0]},)`
      : `(${shape.join(', ')},)`;
    const header = `{'descr': '${dtype}', 'fortran_order': False, 'shape': ${shapeStr}, }`;

    // Pad header so that (magic + 2 version + 2 headerLen + header + \n) is
    // divisible by 64
    const prefix = 10; // 6 magic + 2 version + 2 headerLen
    let padded = header + ' ';
    while ((prefix + padded.length + 1) % 64 !== 0) padded += ' ';
    padded += '\n';

    const headerLen = padded.length;
    const totalHeader = prefix + headerLen;

    const buf = new ArrayBuffer(totalHeader + data.byteLength);
    const view = new DataView(buf);

    // Magic string  \x93NUMPY
    view.setUint8(0, 0x93);
    view.setUint8(1, 0x4E); // N
    view.setUint8(2, 0x55); // U
    view.setUint8(3, 0x4D); // M
    view.setUint8(4, 0x50); // P
    view.setUint8(5, 0x59); // Y

    // Version 1.0
    view.setUint8(6, 1);
    view.setUint8(7, 0);

    // HEADER_LEN (little-endian uint16)
    view.setUint16(8, headerLen, true);

    // Header bytes (ASCII)
    const headerBytes = new Uint8Array(buf, prefix, headerLen);
    for (let i = 0; i < padded.length; i++) {
      headerBytes[i] = padded.charCodeAt(i);
    }

    // Data bytes
    new Uint8Array(buf, totalHeader).set(new Uint8Array(data.buffer, data.byteOffset, data.byteLength));

    return new Blob([buf], { type: 'application/octet-stream' });
  }

  // ─── Public API ───────────────────────────────────────────────────────────

  /**
   * Check whether the MedSAM2 server is reachable.
   * Calls the Gradio /call/health endpoint.
   * Returns { available: boolean, ...serverInfo }
   */
  async function checkHealth() {
    try {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 8000);

      // POST → event_id
      const submitResp = await fetch(`${serverUrl}/call/health`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ data: [] }),
        signal: controller.signal,
      });
      clearTimeout(timeout);

      if (!submitResp.ok) return { available: false };
      const { event_id } = await submitResp.json();

      // GET SSE result
      const data = await readGradioSSE(`${serverUrl}/call/health/${event_id}`);
      const info = JSON.parse(data[0]);
      return { available: true, ...info };
    } catch {
      return { available: false };
    }
  }

  /**
   * Send a 3-D volume + bounding-box prompt to MedSAM2 and receive a
   * binary segmentation mask in return.
   *
   * @param {Float32Array} volumeData  - Flat 3-D volume (x + y*X + z*X*Y)
   * @param {number[]}     shape       - [X, Y, Z]
   * @param {number[]}     spacing     - [sx, sy, sz] in mm
   * @param {Object}       box         - { x1, y1, x2, y2, sliceIdx }
   * @returns {Uint8Array} Binary mask with the same number of elements as the volume
   */
  async function segment(volumeData, shape, spacing, box) {
    const [X, Y, Z] = shape;
    const total = X * Y * Z;

    // ── 1. Window to brain CT range and normalise to 0–255 ──────────────────
    const WL = 40;
    const WW = 80;
    const loHU = WL - WW / 2;  //  0
    const hiHU = WL + WW / 2;  // 80

    const windowed = new Uint8Array(total);
    for (let i = 0; i < total; i++) {
      const hu = volumeData[i];
      const clamped = Math.max(loHU, Math.min(hiHU, hu));
      windowed[i] = Math.round(((clamped - loHU) / WW) * 255);
    }

    // ── 2. Reshape to [Z, Y, X] (server expects [D, H, W]) ─────────────────
    // The web app stores voxels as x + y*X + z*X*Y (XYZ order).
    // NumPy / MedSAM2 expects [depth, height, width] = [Z, Y, X].
    const volumeZYX = new Uint8Array(total);
    for (let z = 0; z < Z; z++) {
      for (let y = 0; y < Y; y++) {
        for (let x = 0; x < X; x++) {
          volumeZYX[z * Y * X + y * X + x] = windowed[x + y * X + z * X * Y];
        }
      }
    }

    // ── 3. Create .npy blob ─────────────────────────────────────────────────
    const npyBlob = createNpyBlob(volumeZYX, [Z, Y, X]);

    // ── 4. Upload .npy to Gradio server ─────────────────────────────────────
    const fileData = await uploadFile(npyBlob, 'volume.npy');

    // ── 5. Build box prompt ─────────────────────────────────────────────────
    const boxObj = {
      x1: box.x1,
      y1: box.y1,
      x2: box.x2,
      y2: box.y2,
      slice_idx: box.sliceIdx !== undefined ? box.sliceIdx : Math.floor(Z / 2),
    };

    // ── 6. Call /segment via Gradio API ─────────────────────────────────────
    const result = await callGradioApi('segment', [
      fileData,
      JSON.stringify(boxObj),
      JSON.stringify([spacing[0], spacing[1], spacing[2]]),
    ]);

    // ── 7. Parse result ─────────────────────────────────────────────────────
    const resultObj = JSON.parse(result[0]);
    if (resultObj.error) throw new Error(resultObj.error);

    // Decode mask: base64 → gzip → npy → Uint8Array
    const maskB64 = resultObj.mask_b64_gzip;
    if (!maskB64) {
      throw new Error('Server response did not include mask_b64_gzip field.');
    }

    const binaryStr = atob(maskB64);
    const maskCompressed = new Uint8Array(binaryStr.length);
    for (let i = 0; i < binaryStr.length; i++) {
      maskCompressed[i] = binaryStr.charCodeAt(i);
    }

    // Decompress gzip
    const decompressed = pako.inflate(maskCompressed);

    // The server wraps the mask in np.save (npy format). We need to skip the
    // npy header to get to the raw uint8 data.
    const maskNpy = new Uint8Array(decompressed);
    const maskData = parseNpyData(maskNpy);

    // ── 8. Transpose mask back from [Z, Y, X] → flat [x + y*X + z*X*Y] ──
    const mask = new Uint8Array(total);
    for (let z = 0; z < Z; z++) {
      for (let y = 0; y < Y; y++) {
        for (let x = 0; x < X; x++) {
          mask[x + y * X + z * X * Y] = maskData[z * Y * X + y * X + x];
        }
      }
    }

    return mask;
  }

  /**
   * Parse raw bytes of a .npy file and return a Uint8Array of the data.
   * Handles the npy header to find the data offset.
   */
  function parseNpyData(bytes) {
    // Magic: \x93NUMPY, then version (2 bytes), then HEADER_LEN (2 or 4 bytes)
    const version = bytes[6];
    let headerLen;
    if (version === 1) {
      headerLen = bytes[8] | (bytes[9] << 8); // uint16 LE
      return new Uint8Array(bytes.buffer, bytes.byteOffset + 10 + headerLen);
    } else if (version === 2) {
      headerLen = bytes[8] | (bytes[9] << 8) | (bytes[10] << 16) | (bytes[11] << 24); // uint32 LE
      return new Uint8Array(bytes.buffer, bytes.byteOffset + 12 + headerLen);
    }
    // Fallback: try to find the newline that ends the header
    let offset = 10;
    while (offset < bytes.length && bytes[offset - 1] !== 0x0A) offset++;
    return new Uint8Array(bytes.buffer, bytes.byteOffset + offset);
  }

  return { setServerUrl, getServerUrl, checkHealth, segment };

})();
