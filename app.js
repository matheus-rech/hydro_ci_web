/**
 * HydroMorph — Hydrocephalus Morphometrics Pipeline
 * Runs entirely in the browser. No data leaves your device.
 * Author: Matheus Machado Rech
 */

'use strict';

// ─── State ──────────────────────────────────────────────────────────────────
let appState = {
  screen: 'upload',        // 'upload' | 'processing' | 'results'
  volume: null,            // { shape, spacing, affine, data, header }
  mask: null,              // Uint8Array ventricle mask
  results: null,           // computed metrics
  currentAxialSlice: 0,
  showMask: true,
  fileName: '',
  fileSize: 0,
  medsam2Available: false, // true when MedSAM2 server is reachable
  medsam2Box: null,        // { x1, y1, x2, y2, sliceIdx } — optional prompt
  segmentationMethod: 'threshold', // 'threshold' | 'medsam2'
};

// ─── Volume Helpers ──────────────────────────────────────────────────────────
function getVoxel(data, shape, x, y, z) {
  if (x < 0 || y < 0 || z < 0 || x >= shape[0] || y >= shape[1] || z >= shape[2]) return 0;
  return data[x + y * shape[0] + z * shape[0] * shape[1]];
}

function setVoxel(data, shape, x, y, z, val) {
  if (x < 0 || y < 0 || z < 0 || x >= shape[0] || y >= shape[1] || z >= shape[2]) return;
  data[x + y * shape[0] + z * shape[0] * shape[1]] = val;
}

function voxelIndex(shape, x, y, z) {
  return x + y * shape[0] + z * shape[0] * shape[1];
}

// ─── Morphological Operations ────────────────────────────────────────────────

/**
 * 3D erosion with 6-connectivity (box kernel 3x3x3 but checking 6 neighbors)
 * For speed, uses 6-connectivity rather than full box for erosion
 */
function erode3D(mask, shape, iterations = 1) {
  const [X, Y, Z] = shape;
  const total = X * Y * Z;
  let src = new Uint8Array(mask);
  let dst = new Uint8Array(total);

  const neighbors = [
    [-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]
  ];

  for (let iter = 0; iter < iterations; iter++) {
    dst.fill(0);
    for (let z = 1; z < Z-1; z++) {
      for (let y = 1; y < Y-1; y++) {
        for (let x = 1; x < X-1; x++) {
          if (src[voxelIndex(shape, x, y, z)] === 0) continue;
          let keep = true;
          for (const [dx,dy,dz] of neighbors) {
            if (src[voxelIndex(shape, x+dx, y+dy, z+dz)] === 0) {
              keep = false; break;
            }
          }
          if (keep) dst[voxelIndex(shape, x, y, z)] = 1;
        }
      }
    }
    src = dst;
    dst = new Uint8Array(total);
  }
  return src;
}

/**
 * 3D dilation with 6-connectivity
 */
function dilate3D(mask, shape, iterations = 1) {
  const [X, Y, Z] = shape;
  const total = X * Y * Z;
  let src = new Uint8Array(mask);
  let dst = new Uint8Array(total);

  const neighbors = [
    [-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]
  ];

  for (let iter = 0; iter < iterations; iter++) {
    dst.set(src);
    for (let z = 1; z < Z-1; z++) {
      for (let y = 1; y < Y-1; y++) {
        for (let x = 1; x < X-1; x++) {
          if (src[voxelIndex(shape, x, y, z)] === 0) continue;
          for (const [dx,dy,dz] of neighbors) {
            dst[voxelIndex(shape, x+dx, y+dy, z+dz)] = 1;
          }
        }
      }
    }
    src = dst;
    dst = new Uint8Array(total);
  }
  return src;
}

// Opening = erode then dilate
function opening3D(mask, shape, iterations = 1) {
  return dilate3D(erode3D(mask, shape, iterations), shape, iterations);
}

// Closing = dilate then erode
function closing3D(mask, shape, iterations = 1) {
  return erode3D(dilate3D(mask, shape, iterations), shape, iterations);
}

// ─── Connected Components (BFS) ──────────────────────────────────────────────

/**
 * 3D connected components with 6-connectivity.
 * Returns { labels: Int32Array, counts: Map<label, count>, numLabels }
 */
function connectedComponents3D(mask, shape) {
  const [X, Y, Z] = shape;
  const total = X * Y * Z;
  const labels = new Int32Array(total);
  const counts = new Map();
  let nextLabel = 1;

  // Pre-allocated BFS queue (flat triplet encoding)
  const queue = new Int32Array(Math.min(total * 3, 6 * 1024 * 1024)); // max 2M voxels queue
  
  const neighbors = [
    [-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]
  ];

  for (let z = 0; z < Z; z++) {
    for (let y = 0; y < Y; y++) {
      for (let x = 0; x < X; x++) {
        const idx = voxelIndex(shape, x, y, z);
        if (mask[idx] === 0 || labels[idx] !== 0) continue;

        const label = nextLabel++;
        labels[idx] = label;
        let count = 1;

        // BFS
        let head = 0, tail = 0;
        queue[tail++] = x;
        queue[tail++] = y;
        queue[tail++] = z;

        while (head < tail) {
          const cx = queue[head++];
          const cy = queue[head++];
          const cz = queue[head++];

          for (const [dx,dy,dz] of neighbors) {
            const nx = cx+dx, ny = cy+dy, nz = cz+dz;
            if (nx < 0 || ny < 0 || nz < 0 || nx >= X || ny >= Y || nz >= Z) continue;
            const nidx = voxelIndex(shape, nx, ny, nz);
            if (mask[nidx] === 0 || labels[nidx] !== 0) continue;
            labels[nidx] = label;
            count++;
            if (tail + 3 < queue.length) {
              queue[tail++] = nx;
              queue[tail++] = ny;
              queue[tail++] = nz;
            }
          }
        }
        counts.set(label, count);
      }
    }
  }
  return { labels, counts, numLabels: nextLabel - 1 };
}

/**
 * Keep only the largest component(s) in a mask.
 */
function keepLargestComponent(mask, shape, minSize = 1) {
  const { labels, counts } = connectedComponents3D(mask, shape);
  const total = mask.length;
  const result = new Uint8Array(total);

  // Find largest
  let maxLabel = -1, maxCount = 0;
  for (const [label, count] of counts) {
    if (count > maxCount) { maxCount = count; maxLabel = label; }
  }

  if (maxLabel === -1 || maxCount < minSize) return result;
  for (let i = 0; i < total; i++) {
    if (labels[i] === maxLabel) result[i] = 1;
  }
  return result;
}

/**
 * Keep all components larger than minSize.
 */
function keepLargeComponents(mask, shape, minSize = 500) {
  const { labels, counts } = connectedComponents3D(mask, shape);
  const total = mask.length;
  const result = new Uint8Array(total);

  const keepLabels = new Set();
  for (const [label, count] of counts) {
    if (count >= minSize) keepLabels.add(label);
  }

  for (let i = 0; i < total; i++) {
    if (keepLabels.has(labels[i])) result[i] = 1;
  }
  return result;
}

// ─── Progress Management ──────────────────────────────────────────────────────

let progressSteps = [];
let currentStep = 0;

function initProgress(steps) {
  progressSteps = steps;
  currentStep = 0;
  renderProgressSteps();
}

function advanceProgress(stepIndex, message) {
  currentStep = stepIndex;
  renderProgressSteps();
  if (message) {
    const el = document.getElementById('progress-detail');
    if (el) el.textContent = message;
  }
}

function renderProgressSteps() {
  const container = document.getElementById('progress-steps');
  if (!container) return;
  container.innerHTML = progressSteps.map((step, i) => {
    let cls = 'step-pending';
    let icon = '○';
    if (i < currentStep) { cls = 'step-done'; icon = '✓'; }
    else if (i === currentStep) { cls = 'step-active'; icon = '◉'; }
    return `<div class="progress-step ${cls}"><span class="step-icon">${icon}</span><span class="step-label">${step}</span></div>`;
  }).join('');
}

// ─── Main Pipeline ────────────────────────────────────────────────────────────

async function delay(ms = 0) {
  return new Promise(r => setTimeout(r, ms));
}

async function runPipeline(volume) {
  const steps = [
    'Parsing scan header',
    'Building brain mask',
    'Extracting CSF voxels',
    'Morphological filtering',
    'Isolating ventricles',
    'Computing Evans Index',
    'Computing callosal angle',
    'Computing volume',
    'Generating report'
  ];
  initProgress(steps);
  await delay(30);

  const { shape, spacing, data } = volume;
  const [X, Y, Z] = shape;
  const total = X * Y * Z;

  // ── Step 0: Header already parsed ──────────────────────────────────────────
  advanceProgress(0, `Volume: ${X}×${Y}×${Z}, spacing: ${spacing.map(s=>s.toFixed(2)).join('×')} mm`);
  updateVolumeMetadata(volume);
  await delay(50);

  // Check MedSAM2 server availability in the background (non-blocking)
  MedSAMClient.checkHealth().then(result => {
    appState.medsam2Available = result.available;
    updateServerStatusUI(result.available ? 'connected' : null);
  }).catch(() => {
    appState.medsam2Available = false;
  });

  // ── Step 1: Brain mask ─────────────────────────────────────────────────────
  advanceProgress(1, 'Thresholding brain tissue (HU: -5 to 80)...');
  await delay(10);

  const brainMaskRaw = new Uint8Array(total);
  for (let i = 0; i < total; i++) {
    const hu = data[i];
    brainMaskRaw[i] = (hu >= -5 && hu <= 80) ? 1 : 0;
  }

  // Closing to fill small gaps (2 iterations)
  advanceProgress(1, 'Closing brain mask...');
  await delay(10);
  let brainMask = closing3D(brainMaskRaw, shape, 2);

  // Keep largest component
  advanceProgress(1, 'Keeping largest brain component...');
  await delay(10);
  brainMask = keepLargestComponent(brainMask, shape, 1000);

  const brainVoxCount = brainMask.reduce((a, v) => a + v, 0);
  advanceProgress(1, `Brain mask: ${brainVoxCount.toLocaleString()} voxels`);
  await delay(20);

  // ── Step 2: CSF mask ───────────────────────────────────────────────────────
  advanceProgress(2, 'Extracting CSF (HU: 0 to 22) within brain...');
  await delay(10);

  const csfMask = new Uint8Array(total);
  for (let i = 0; i < total; i++) {
    const hu = data[i];
    csfMask[i] = (brainMask[i] === 1 && hu >= 0 && hu <= 22) ? 1 : 0;
  }

  const csfCount = csfMask.reduce((a, v) => a + v, 0);
  advanceProgress(2, `CSF voxels: ${csfCount.toLocaleString()}`);
  await delay(20);

  // ── Step 3: Morphological filtering (denoise) ────────────────────────────
  advanceProgress(3, 'Applying morphological filtering...');
  await delay(10);
  // For high-res volumes (small spacing), erosion removes too much.
  // For low-res volumes (large spacing), erosion is too coarse.
  // Only apply opening when spacing is in a reasonable range (0.7-2.5mm).
  const minSpacingXY = Math.min(spacing[0], spacing[1]);
  let ventMask;
  if (minSpacingXY < 0.7 || minSpacingXY > 2.5) {
    // Skip morphological opening, rely on component filtering
    advanceProgress(3, 'Adaptive filtering — skipping erosion for this resolution...');
    ventMask = new Uint8Array(csfMask);
  } else {
    ventMask = opening3D(csfMask, shape, 1);
  }
  await delay(20);

  // ── Step 4: Restrict to central 60% ───────────────────────────────────────
  advanceProgress(4, 'Restricting to central brain region (ventricles are central)...');
  await delay(10);

  // Compute brain bounding box
  let bxMin=X, bxMax=0, byMin=Y, byMax=0, bzMin=Z, bzMax=0;
  for (let z = 0; z < Z; z++) for (let y = 0; y < Y; y++) for (let x = 0; x < X; x++) {
    if (brainMask[voxelIndex(shape, x, y, z)] === 1) {
      if (x < bxMin) bxMin = x; if (x > bxMax) bxMax = x;
      if (y < byMin) byMin = y; if (y > byMax) byMax = y;
      if (z < bzMin) bzMin = z; if (z > bzMax) bzMax = z;
    }
  }

  // Central 60%: trim 20% from each side
  const marginX = Math.floor((bxMax - bxMin) * 0.20);
  const marginY = Math.floor((byMax - byMin) * 0.20);
  const marginZ = Math.floor((bzMax - bzMin) * 0.10); // less z margin — ventricles span more vertically

  const cropXmin = bxMin + marginX, cropXmax = bxMax - marginX;
  const cropYmin = byMin + marginY, cropYmax = byMax - marginY;
  const cropZmin = bzMin + marginZ, cropZmax = bzMax - marginZ;

  // Apply central crop restriction
  for (let z = 0; z < Z; z++) for (let y = 0; y < Y; y++) for (let x = 0; x < X; x++) {
    const idx = voxelIndex(shape, x, y, z);
    if (ventMask[idx] === 0) continue;
    if (x < cropXmin || x > cropXmax || y < cropYmin || y > cropYmax ||
        z < cropZmin || z > cropZmax) {
      ventMask[idx] = 0;
    }
  }

  // ── Step 4b: Keep components > threshold (adaptive to resolution) ─────────
  // Scale minimum component size with voxel volume so threshold is ~1 mL
  const voxVol = spacing[0] * spacing[1] * spacing[2];
  const minVolumeMl = 0.5; // minimum 0.5 mL
  const minComponentSize = Math.max(50, Math.round((minVolumeMl * 1000) / voxVol));
  advanceProgress(4, `Filtering connected components (>${minComponentSize} voxels, ~${minVolumeMl} mL)...`);
  await delay(10);
  ventMask = keepLargeComponents(ventMask, shape, minComponentSize);

  // ── Optional: MedSAM2 AI segmentation override ─────────────────────────────
  if (appState.medsam2Available && appState.medsam2Box) {
    advanceProgress(4, 'Running MedSAM2 AI segmentation...');
    try {
      const aiMask = await MedSAMClient.segment(data, shape, spacing, appState.medsam2Box);
      // Validate returned mask size
      if (aiMask && aiMask.length === total) {
        ventMask = aiMask;
        appState.segmentationMethod = 'medsam2';
        advanceProgress(4, 'MedSAM2 segmentation applied.');
      } else {
        throw new Error('Mask size mismatch from server.');
      }
    } catch (err) {
      advanceProgress(4, `MedSAM2 failed (${err.message}) — using threshold fallback.`);
      appState.segmentationMethod = 'threshold';
    }
  } else {
    appState.segmentationMethod = 'threshold';
  }

  const ventCount = ventMask.reduce((a, v) => a + v, 0);
  advanceProgress(4, `Ventricle voxels: ${ventCount.toLocaleString()} [${appState.segmentationMethod}]`);
  await delay(20);

  if (ventCount < 100) {
    throw new Error(`Very few ventricle voxels found (${ventCount}). Is this a head CT with HU values? Check that your file is a CT scan in Hounsfield Units.`);
  }

  // ── Step 5: Evans Index ────────────────────────────────────────────────────
  advanceProgress(5, 'Computing Evans Index per axial slice...');
  await delay(10);

  const evansResult = computeEvansIndex(data, ventMask, shape, spacing);
  await delay(20);

  // ── Step 6: Callosal Angle ─────────────────────────────────────────────────
  advanceProgress(6, 'Computing callosal angle on coronal view...');
  await delay(10);

  const callosalResult = computeCallosalAngle(ventMask, shape, spacing);
  await delay(20);

  // ── Step 7: Volume ─────────────────────────────────────────────────────────
  advanceProgress(7, 'Computing ventricle volume...');
  const voxelVol = spacing[0] * spacing[1] * spacing[2]; // mm³
  const ventVolMm3 = ventCount * voxelVol;
  const ventVolMl = ventVolMm3 / 1000;
  await delay(10);

  // ── Step 8: NPH Probability ────────────────────────────────────────────────
  advanceProgress(8, 'Generating clinical report...');
  await delay(20);

  let nphScore = 0;
  if (evansResult.maxEvans > 0.3) nphScore++;
  if (callosalResult.angleDeg !== null && callosalResult.angleDeg < 90) nphScore++;
  if (ventVolMl > 50) nphScore++;
  const nphPct = Math.round((nphScore / 3) * 100);

  const results = {
    evansIndex: evansResult.maxEvans,
    evansSlice: evansResult.bestSlice,
    evansData: evansResult,
    callosalAngle: callosalResult.angleDeg,
    callosalSlice: callosalResult.bestCoronalSlice,
    callosalData: callosalResult,
    ventVolMl,
    ventVolMm3,
    nphScore,
    nphPct,
    brainBbox: { bxMin, bxMax, byMin, byMax, bzMin, bzMax },
    ventCount,
    brainVoxCount,
    shape,
    spacing
  };

  appState.mask = ventMask;
  appState.results = results;
  appState.currentAxialSlice = evansResult.bestSlice >= 0 ? evansResult.bestSlice : Math.floor(Z/2);

  return results;
}

// ─── Evans Index Computation ──────────────────────────────────────────────────

function computeEvansIndex(data, ventMask, shape, spacing) {
  const [X, Y, Z] = shape;
  let maxEvans = 0;
  let bestSlice = -1;
  const perSlice = [];

  for (let z = 0; z < Z; z++) {
    // Find all ventricle voxels on this slice
    let ventLeft = X, ventRight = 0, ventTop = 0, ventCount = 0;
    for (let y = 0; y < Y; y++) {
      for (let x = 0; x < X; x++) {
        if (ventMask[voxelIndex(shape, x, y, z)] === 1) {
          ventCount++;
          if (x < ventLeft) ventLeft = x;
          if (x > ventRight) ventRight = x;
          if (y > ventTop) ventTop = y;
        }
      }
    }

    if (ventCount < 20) continue;

    const ventWidthMm = (ventRight - ventLeft) * spacing[0];

    // Find skull width: use bone HU > 300 or fallback to brain extent
    let skullLeft = X, skullRight = 0;

    // Try bone thresholding first
    let boneCount = 0;
    for (let y = 0; y < Y; y++) {
      for (let x = 0; x < X; x++) {
        const hu = getVoxel(data, shape, x, y, z);
        if (hu > 300) {
          boneCount++;
          if (x < skullLeft) skullLeft = x;
          if (x > skullRight) skullRight = x;
        }
      }
    }

    // Fallback to soft tissue extent if few bone voxels
    if (boneCount < 10 || (skullRight - skullLeft) < 50) {
      skullLeft = X; skullRight = 0;
      for (let y = 0; y < Y; y++) {
        for (let x = 0; x < X; x++) {
          const hu = getVoxel(data, shape, x, y, z);
          if (hu > -20 && hu < 1000) {
            if (x < skullLeft) skullLeft = x;
            if (x > skullRight) skullRight = x;
          }
        }
      }
    }

    if (skullRight <= skullLeft) continue;
    const skullWidthMm = (skullRight - skullLeft) * spacing[0];

    if (skullWidthMm < 50) continue; // sanity check

    const evans = ventWidthMm / skullWidthMm;

    perSlice.push({ z, evans, ventWidthMm, skullWidthMm, ventLeft, ventRight, skullLeft, skullRight });

    if (evans > maxEvans) {
      maxEvans = evans;
      bestSlice = z;
    }
  }

  return { maxEvans, bestSlice, perSlice };
}

// ─── Callosal Angle Computation ───────────────────────────────────────────────

function computeCallosalAngle(ventMask, shape, spacing) {
  const [X, Y, Z] = shape;

  // Find coronal slice (axis 1 = Y) with largest ventricle cross-section
  let maxCount = 0, bestY = -1;
  for (let y = 0; y < Y; y++) {
    let count = 0;
    for (let z = 0; z < Z; z++) for (let x = 0; x < X; x++) {
      if (ventMask[voxelIndex(shape, x, y, z)] === 1) count++;
    }
    if (count > maxCount) { maxCount = count; bestY = y; }
  }

  if (bestY === -1 || maxCount < 20) {
    return { angleDeg: null, bestCoronalSlice: -1, vertex: null, leftPt: null, rightPt: null };
  }

  // On the coronal slice (x vs z), find the top (highest z = superior)
  let topZ = 0;
  let leftX = X, rightX = 0;
  const midX = Math.floor(X / 2);

  // Find topmost voxels on left and right halves
  let topZLeft = 0, topZRight = 0;
  let leftPtX = -1, leftPtZ = -1, rightPtX = -1, rightPtZ = -1;

  for (let z = 0; z < Z; z++) {
    for (let x = 0; x < X; x++) {
      if (ventMask[voxelIndex(shape, x, bestY, z)] === 1) {
        if (z > topZ) topZ = z;
        if (x < leftX) leftX = x;
        if (x > rightX) rightX = x;
        if (x < midX && z > topZLeft) { topZLeft = z; leftPtX = x; leftPtZ = z; }
        if (x >= midX && z > topZRight) { topZRight = z; rightPtX = x; rightPtZ = z; }
      }
    }
  }

  // Vertex = topmost point overall (centroid of top voxels)
  let vertexSumX = 0, vertexSumZ = 0, vertexN = 0;
  const topZThresh = topZ - 3;
  for (let z = topZThresh; z <= topZ; z++) {
    for (let x = 0; x < X; x++) {
      if (ventMask[voxelIndex(shape, x, bestY, z)] === 1) {
        vertexSumX += x; vertexSumZ += z; vertexN++;
      }
    }
  }
  if (vertexN === 0) {
    return { angleDeg: null, bestCoronalSlice: bestY, vertex: null, leftPt: null, rightPt: null };
  }
  const vx = vertexSumX / vertexN;
  const vz = vertexSumZ / vertexN;

  // Left and right points: bottom-most extremes on each side
  let bLeftX=-1, bLeftZ=Z, bRightX=-1, bRightZ=Z;
  for (let z = 0; z < Z; z++) {
    for (let x = 0; x < X; x++) {
      if (ventMask[voxelIndex(shape, x, bestY, z)] === 1) {
        if (x < midX && z < bLeftZ) { bLeftZ = z; bLeftX = x; }
        if (x >= midX && z < bRightZ) { bRightZ = z; bRightX = x; }
      }
    }
  }

  if (bLeftX < 0 || bRightX < 0) {
    return { angleDeg: null, bestCoronalSlice: bestY, vertex: null, leftPt: null, rightPt: null };
  }

  // Vectors from vertex to bottom-left and bottom-right (correcting for spacing)
  const lx = (bLeftX - vx) * spacing[0];
  const lz = (bLeftZ - vz) * spacing[2];
  const rx = (bRightX - vx) * spacing[0];
  const rz = (bRightZ - vz) * spacing[2];

  const dotProd = lx*rx + lz*rz;
  const magL = Math.sqrt(lx*lx + lz*lz);
  const magR = Math.sqrt(rx*rx + rz*rz);

  if (magL < 0.001 || magR < 0.001) {
    return { angleDeg: null, bestCoronalSlice: bestY, vertex: null, leftPt: null, rightPt: null };
  }

  const cosAngle = Math.max(-1, Math.min(1, dotProd / (magL * magR)));
  const angleDeg = Math.round(Math.acos(cosAngle) * (180 / Math.PI));

  return {
    angleDeg,
    bestCoronalSlice: bestY,
    vertex: { x: vx, z: vz },
    leftPt: { x: bLeftX, z: bLeftZ },
    rightPt: { x: bRightX, z: bRightZ },
    midX
  };
}

// ─── Canvas Rendering ─────────────────────────────────────────────────────────

function renderAxialSlice(canvas, volume, mask, sliceIndex, showMask) {
  const { shape, data } = volume;
  const [X, Y] = shape;
  const ctx = canvas.getContext('2d');

  // Scale canvas to fit container while maintaining aspect ratio
  const containerW = canvas.parentElement ? canvas.parentElement.clientWidth : 320;
  const scale = Math.min(containerW / X, containerW / Y, 2.0);
  canvas.width = Math.round(X * scale);
  canvas.height = Math.round(Y * scale);
  canvas.style.width = canvas.width + 'px';
  canvas.style.height = canvas.height + 'px';

  const offscreenCanvas = document.createElement('canvas');
  offscreenCanvas.width = X;
  offscreenCanvas.height = Y;
  const offCtx = offscreenCanvas.getContext('2d');
  const imageData = offCtx.createImageData(X, Y);

  const lo = 0, hi = 80;

  for (let y = 0; y < Y; y++) {
    for (let x = 0; x < X; x++) {
      const hu = getVoxel(data, shape, x, y, sliceIndex);
      let gray = Math.floor(((Math.min(Math.max(hu, lo), hi) - lo) / (hi - lo)) * 255);
      const pixIdx = (y * X + x) * 4;

      const isMask = showMask && mask && mask[voxelIndex(shape, x, y, sliceIndex)] === 1;
      if (isMask) {
        imageData.data[pixIdx]   = Math.floor(gray * 0.4 + 88 * 0.6);
        imageData.data[pixIdx+1] = Math.floor(gray * 0.4 + 166 * 0.6);
        imageData.data[pixIdx+2] = Math.floor(gray * 0.4 + 255 * 0.6);
      } else {
        imageData.data[pixIdx]   = gray;
        imageData.data[pixIdx+1] = gray;
        imageData.data[pixIdx+2] = gray;
      }
      imageData.data[pixIdx+3] = 255;
    }
  }

  offCtx.putImageData(imageData, 0, 0);
  ctx.save();
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offscreenCanvas, 0, 0, canvas.width, canvas.height);
  ctx.restore();

  // Draw Evans index annotation on best slice
  if (appState.results && sliceIndex === appState.results.evansSlice) {
    drawEvansAnnotation(ctx, canvas.width, canvas.height, scale);
  }
}

function drawEvansAnnotation(ctx, W, H, scale) {
  const results = appState.results;
  if (!results || !results.evansData) return;
  const perSlice = results.evansData.perSlice;
  const sliceData = perSlice.find(s => s.z === results.evansSlice);
  if (!sliceData) return;

  const { ventLeft, ventRight, skullLeft, skullRight } = sliceData;
  const { shape, spacing } = appState.volume;
  const [X, Y] = shape;
  const sy = Math.floor(Y * 0.4) * scale; // draw at ~40% height

  ctx.lineWidth = 2;
  ctx.strokeStyle = '#58a6ff';
  ctx.beginPath();
  ctx.moveTo(ventLeft * scale, sy);
  ctx.lineTo(ventRight * scale, sy);
  ctx.stroke();

  ctx.strokeStyle = '#ff6e40';
  ctx.beginPath();
  ctx.moveTo(skullLeft * scale, sy + 8);
  ctx.lineTo(skullRight * scale, sy + 8);
  ctx.stroke();

  // Labels
  ctx.fillStyle = '#58a6ff';
  ctx.font = `${Math.max(10, 12*scale)}px JetBrains Mono, monospace`;
  ctx.fillText('V', ventLeft * scale, sy - 4);
  ctx.fillStyle = '#ff6e40';
  ctx.fillText('S', skullLeft * scale, sy + 22);
}

function renderCoronalSlice(canvas, volume, mask, sliceY, callosalData) {
  const { shape, data } = volume;
  const [X, Y, Z] = shape;
  const ctx = canvas.getContext('2d');

  const containerW = canvas.parentElement ? canvas.parentElement.clientWidth : 320;
  const scale = Math.min(containerW / X, containerW / Z, 2.0);
  canvas.width = Math.round(X * scale);
  canvas.height = Math.round(Z * scale);
  canvas.style.width = canvas.width + 'px';
  canvas.style.height = canvas.height + 'px';

  const offscreenCanvas = document.createElement('canvas');
  offscreenCanvas.width = X;
  offscreenCanvas.height = Z;
  const offCtx = offscreenCanvas.getContext('2d');
  const imageData = offCtx.createImageData(X, Z);

  const lo = 0, hi = 80;

  for (let z = 0; z < Z; z++) {
    for (let x = 0; x < X; x++) {
      const hu = getVoxel(data, shape, x, sliceY, z);
      let gray = Math.floor(((Math.min(Math.max(hu, lo), hi) - lo) / (hi - lo)) * 255);
      // Flip z so inferior is at bottom (z=0 -> bottom, z=Z-1 -> top)
      const dispZ = Z - 1 - z;
      const pixIdx = (dispZ * X + x) * 4;
      const isMask = mask && mask[voxelIndex(shape, x, sliceY, z)] === 1;
      if (isMask) {
        imageData.data[pixIdx]   = Math.floor(gray * 0.4 + 88 * 0.6);
        imageData.data[pixIdx+1] = Math.floor(gray * 0.4 + 166 * 0.6);
        imageData.data[pixIdx+2] = Math.floor(gray * 0.4 + 255 * 0.6);
      } else {
        imageData.data[pixIdx]   = gray;
        imageData.data[pixIdx+1] = gray;
        imageData.data[pixIdx+2] = gray;
      }
      imageData.data[pixIdx+3] = 255;
    }
  }

  offCtx.putImageData(imageData, 0, 0);
  ctx.save();
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offscreenCanvas, 0, 0, canvas.width, canvas.height);
  ctx.restore();

  // Draw callosal angle annotation
  if (callosalData && callosalData.vertex && callosalData.leftPt && callosalData.rightPt) {
    drawCallosalAnnotation(ctx, callosalData, Z, scale);
  }
}

function drawCallosalAnnotation(ctx, callosalData, Z, scale) {
  const { vertex, leftPt, rightPt } = callosalData;
  // Flip z
  const vx = vertex.x * scale;
  const vz = (Z - 1 - vertex.z) * scale;
  const lx = leftPt.x * scale;
  const lz = (Z - 1 - leftPt.z) * scale;
  const rx = rightPt.x * scale;
  const rz = (Z - 1 - rightPt.z) * scale;

  // Draw vertex dot
  ctx.fillStyle = '#00d4d4';
  ctx.beginPath();
  ctx.arc(vx, vz, 5, 0, Math.PI*2);
  ctx.fill();

  // Draw vectors
  ctx.strokeStyle = '#00d4d4';
  ctx.lineWidth = 2.5;
  ctx.setLineDash([4, 3]);
  ctx.beginPath();
  ctx.moveTo(vx, vz); ctx.lineTo(lx, lz); ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(vx, vz); ctx.lineTo(rx, rz); ctx.stroke();
  ctx.setLineDash([]);

  // Angle arc
  const angle = callosalData.angleDeg;
  if (angle !== null) {
    ctx.fillStyle = '#00d4d4';
    ctx.font = `bold 14px JetBrains Mono, monospace`;
    ctx.fillText(`${angle}°`, vx + 8, vz - 8);
  }

  // Left/right dots
  ctx.fillStyle = '#ff6e40';
  ctx.beginPath(); ctx.arc(lx, lz, 4, 0, Math.PI*2); ctx.fill();
  ctx.beginPath(); ctx.arc(rx, rz, 4, 0, Math.PI*2); ctx.fill();
}

// ─── UI Metadata Update ───────────────────────────────────────────────────────

function updateVolumeMetadata(volume) {
  const el = document.getElementById('volume-metadata');
  if (!el) return;
  const { shape, spacing, header } = volume;
  const fileSizeMB = (appState.fileSize / 1024 / 1024).toFixed(1);
  el.innerHTML = `
    <div class="meta-item"><span class="meta-label">Shape</span><span class="meta-value">${shape[0]}×${shape[1]}×${shape[2]}</span></div>
    <div class="meta-item"><span class="meta-label">Spacing</span><span class="meta-value">${spacing[0].toFixed(2)}×${spacing[1].toFixed(2)}×${spacing[2].toFixed(2)} mm</span></div>
    <div class="meta-item"><span class="meta-label">Datatype</span><span class="meta-value">INT${header.bitpix}</span></div>
    <div class="meta-item"><span class="meta-label">File size</span><span class="meta-value">${fileSizeMB} MB</span></div>
  `;
}

// ─── Screen Transitions ───────────────────────────────────────────────────────

function showScreen(name) {
  document.querySelectorAll('.screen').forEach(s => s.classList.remove('active'));
  const el = document.getElementById(`screen-${name}`);
  if (el) {
    el.classList.add('active');
    appState.screen = name;
  }
}

// ─── File Handling ────────────────────────────────────────────────────────────

/**
 * Main entry point: route one or more files to the correct handler.
 */
async function handleFiles(files) {
  if (!files || files.length === 0) return;

  const fileArray = Array.from(files);
  const firstFile = fileArray[0];
  const name = firstFile.name.toLowerCase();

  // NIfTI
  if (name.endsWith('.nii') || name.endsWith('.nii.gz')) {
    return handleNiftiFile(firstFile);
  }

  // Explicit DICOM extension
  if (name.endsWith('.dcm') || firstFile.type === 'application/dicom') {
    return handleDicomFiles(fileArray);
  }

  // Multiple files with no clear NIfTI extension → assume DICOM series
  if (fileArray.length > 1) {
    return handleDicomFiles(fileArray);
  }

  // 2D image formats
  if (name.endsWith('.png') || name.endsWith('.jpg') || name.endsWith('.jpeg') || name.endsWith('.bmp')) {
    return handleImageFile(firstFile);
  }

  // Try to detect DICOM by magic bytes ("DICM" at offset 128)
  try {
    const header = new Uint8Array(await firstFile.slice(0, 132).arrayBuffer());
    if (header.length >= 132 &&
        String.fromCharCode(header[128], header[129], header[130], header[131]) === 'DICM') {
      return handleDicomFiles([firstFile]);
    }
  } catch { /* ignore */ }

  showError('Unsupported file format. Please use NIfTI (.nii/.nii.gz), DICOM (.dcm), or image files (.png/.jpg).');
}

/**
 * Handle a NIfTI file (.nii / .nii.gz).
 */
async function handleNiftiFile(file) {
  appState.fileName = file.name;
  appState.fileSize = file.size;
  showScreen('processing');

  const steps = [
    'Parsing NIfTI header', 'Building brain mask', 'Extracting CSF voxels',
    'Morphological filtering', 'Isolating ventricles', 'Computing Evans Index',
    'Computing callosal angle', 'Computing volume', 'Generating report'
  ];
  initProgress(steps);

  const fnEl = document.getElementById('processing-filename');
  if (fnEl) fnEl.textContent = file.name;

  try {
    advanceProgress(0, 'Reading file...');
    const arrayBuffer = await file.arrayBuffer();

    advanceProgress(0, 'Decompressing & parsing NIfTI...');
    await delay(30);
    const volume = await NiftiReader.parse(arrayBuffer);
    appState.volume = volume;

    const results = await runPipeline(volume);
    buildResultsUI(results);
    showScreen('results');
  } catch (err) {
    console.error('NIfTI pipeline error:', err);
    showError(err.message || 'An error occurred during processing.');
  }
}

/**
 * Handle one or more DICOM files (.dcm series).
 */
async function handleDicomFiles(files) {
  appState.fileName = files.length > 1
    ? `DICOM series (${files.length} files)`
    : files[0].name;
  appState.fileSize = files.reduce((sum, f) => sum + f.size, 0);
  showScreen('processing');

  const fnEl = document.getElementById('processing-filename');
  if (fnEl) fnEl.textContent = appState.fileName;

  try {
    advanceProgress(0, `Reading ${files.length} DICOM file(s)...`);
    const arrayBuffers = await Promise.all(files.map(f => f.arrayBuffer()));

    advanceProgress(0, 'Parsing DICOM headers and pixel data...');
    await delay(30);
    const volume = await DicomReader.parseSeries(arrayBuffers);
    appState.volume = volume;

    const results = await runPipeline(volume);
    buildResultsUI(results);
    showScreen('results');
  } catch (err) {
    console.error('DICOM pipeline error:', err);
    showError(err.message || 'Failed to parse DICOM files.');
  }
}

/**
 * Handle a single 2D image file (PNG/JPG/BMP).
 */
async function handleImageFile(file) {
  appState.fileName = file.name;
  appState.fileSize = file.size;
  showScreen('processing');

  const fnEl = document.getElementById('processing-filename');
  if (fnEl) fnEl.textContent = file.name;

  try {
    advanceProgress(0, 'Loading image...');
    const volume = await DicomReader.parseImage(file);
    appState.volume = volume;

    advanceProgress(0, 'Single 2D image — running limited analysis. Connect MedSAM2 for AI segmentation.');
    const results = await runPipeline(volume);
    buildResultsUI(results);
    showScreen('results');
  } catch (err) {
    console.error('Image pipeline error:', err);
    showError(err.message || 'Failed to process image.');
  }
}

// Backward-compatible alias so any legacy callers still work
function handleFile(file) {
  return handleFiles([file]);
}


function showError(msg) {
  showScreen('upload');
  const errEl = document.getElementById('upload-error');
  if (errEl) {
    errEl.textContent = msg;
    errEl.style.display = 'block';
    setTimeout(() => { errEl.style.display = 'none'; }, 8000);
  }
}

// ─── Results UI Builder ───────────────────────────────────────────────────────

function buildResultsUI(results) {
  const { evansIndex, callosalAngle, ventVolMl, nphScore, nphPct } = results;

  // NPH Badge
  const badgeEl = document.getElementById('nph-badge');
  if (badgeEl) {
    let label, cls;
    if (nphScore >= 2) { label = 'HIGH'; cls = 'badge-high'; }
    else if (nphScore === 1) { label = 'MODERATE'; cls = 'badge-moderate'; }
    else { label = 'LOW'; cls = 'badge-low'; }
    badgeEl.className = `nph-badge ${cls}`;
    badgeEl.innerHTML = `
      <div class="badge-label">${label}</div>
      <div class="badge-sub">NPH Probability</div>
      <div class="badge-pct">${nphPct}%</div>
      <div class="badge-score">${nphScore}/3 criteria met</div>
    `;
  }

  // Metrics Cards
  setMetricCard('card-evans', evansIndex !== undefined ? evansIndex.toFixed(3) : '—',
    'Evans Index', '>0.3 = abnormal', evansIndex > 0.3 ? 'abnormal' : 'normal');
  setMetricCard('card-angle', callosalAngle !== null ? `${callosalAngle}°` : '—',
    'Callosal Angle', '<90° = abnormal', callosalAngle !== null && callosalAngle < 90 ? 'abnormal' : 'normal');
  setMetricCard('card-volume', `${ventVolMl.toFixed(1)} mL`,
    'Ventricle Volume', '>50 mL = abnormal', ventVolMl > 50 ? 'abnormal' : 'normal');
  setMetricCard('card-nph', `${nphPct}%`,
    'NPH Probability', `${nphScore}/3 criteria`, nphScore >= 2 ? 'abnormal' : (nphScore === 1 ? 'moderate' : 'normal'));

  // Segmentation method indicator
  const methodEl = document.getElementById('segmentation-method-badge');
  if (methodEl) {
    const isMedSAM = appState.segmentationMethod === 'medsam2';
    methodEl.textContent = isMedSAM ? 'AI (MedSAM2)' : 'Threshold';
    methodEl.className = 'method-badge ' + (isMedSAM ? 'method-ai' : 'method-threshold');
    methodEl.style.display = 'inline-flex';
  }

  // Detailed measurements table
  buildMeasurementsTable(results);

  // Sanity checks
  buildSanityChecks(results);

  // Slice viewer setup
  setupSliceViewers(results);
}

function setMetricCard(id, value, label, ref, status) {
  const el = document.getElementById(id);
  if (!el) return;
  el.className = `metric-card metric-${status}`;
  el.innerHTML = `
    <div class="metric-value mono">${value}</div>
    <div class="metric-label">${label}</div>
    <div class="metric-ref">${ref}</div>
  `;
}

function buildMeasurementsTable(results) {
  const el = document.getElementById('measurements-table');
  if (!el) return;
  const { shape, spacing } = results;
  const voxVolMm3 = spacing[0] * spacing[1] * spacing[2];

  el.innerHTML = `
    <table class="data-table">
      <thead><tr><th>Measurement</th><th>Value</th><th>Unit</th><th>Status</th></tr></thead>
      <tbody>
        <tr><td>Evans Index</td><td class="mono">${results.evansIndex.toFixed(4)}</td><td>ratio</td><td class="${results.evansIndex>0.3?'status-abnormal':'status-normal'}">${results.evansIndex>0.3?'ABNORMAL':'NORMAL'}</td></tr>
        <tr><td>Best Evans Slice (axial)</td><td class="mono">${results.evansSlice}</td><td>voxel</td><td>—</td></tr>
        <tr><td>Callosal Angle</td><td class="mono">${results.callosalAngle !== null ? results.callosalAngle : 'N/A'}</td><td>degrees</td><td class="${results.callosalAngle !== null && results.callosalAngle < 90 ? 'status-abnormal':'status-normal'}">${results.callosalAngle !== null ? (results.callosalAngle < 90 ? 'ABNORMAL' : 'NORMAL') : 'N/A'}</td></tr>
        <tr><td>Callosal Slice (coronal)</td><td class="mono">${results.callosalSlice}</td><td>voxel</td><td>—</td></tr>
        <tr><td>Ventricle Volume</td><td class="mono">${results.ventVolMm3.toFixed(0)}</td><td>mm³</td><td>—</td></tr>
        <tr><td>Ventricle Volume</td><td class="mono">${results.ventVolMl.toFixed(2)}</td><td>mL</td><td class="${results.ventVolMl>50?'status-abnormal':'status-normal'}">${results.ventVolMl>50?'ABNORMAL':'NORMAL'}</td></tr>
        <tr><td>Ventricle Voxels</td><td class="mono">${results.ventCount.toLocaleString()}</td><td>voxels</td><td>—</td></tr>
        <tr><td>Voxel Volume</td><td class="mono">${voxVolMm3.toFixed(4)}</td><td>mm³</td><td>—</td></tr>
        <tr><td>Volume (X×Y×Z)</td><td class="mono">${shape[0]}×${shape[1]}×${shape[2]}</td><td>voxels</td><td>—</td></tr>
        <tr><td>Spacing (X×Y×Z)</td><td class="mono">${spacing[0].toFixed(3)}×${spacing[1].toFixed(3)}×${spacing[2].toFixed(3)}</td><td>mm/voxel</td><td>—</td></tr>
        <tr><td>Segmentation Method</td><td class="mono">${appState.segmentationMethod}</td><td>—</td><td>—</td></tr>
      </tbody>
    </table>
  `;
}

function buildSanityChecks(results) {
  const el = document.getElementById('sanity-checks');
  if (!el) return;
  const warnings = [];
  const { evansIndex, callosalAngle, ventVolMl, spacing, shape } = results;

  if (evansIndex > 0.7) warnings.push({ level: 'warn', msg: `Evans Index ${evansIndex.toFixed(3)} is very high (>0.7). Please verify segmentation.` });
  if (evansIndex < 0.1) warnings.push({ level: 'warn', msg: `Evans Index ${evansIndex.toFixed(3)} is very low. Verify ventricles were detected.` });
  if (callosalAngle !== null && callosalAngle > 160) warnings.push({ level: 'warn', msg: `Callosal angle ${callosalAngle}° seems very wide. Verify coronal segmentation.` });
  if (ventVolMl < 5) warnings.push({ level: 'warn', msg: `Ventricle volume ${ventVolMl.toFixed(1)} mL seems very low. Check segmentation.` });
  if (ventVolMl > 200) warnings.push({ level: 'warn', msg: `Ventricle volume ${ventVolMl.toFixed(1)} mL seems very high. Verify segmentation.` });
  if (spacing[0] > 5 || spacing[1] > 5 || spacing[2] > 5) warnings.push({ level: 'warn', msg: `Large voxel spacing detected (${spacing.map(s=>s.toFixed(1)).join('×')} mm). Results may be less accurate.` });

  if (warnings.length === 0) {
    el.innerHTML = `<div class="sanity-ok">✓ All measurements within expected ranges</div>`;
  } else {
    el.innerHTML = warnings.map(w => `<div class="sanity-warn">⚠ ${w.msg}</div>`).join('');
  }
}

function setupSliceViewers(results) {
  const axialSlider = document.getElementById('axial-slider');
  const axialLabel = document.getElementById('axial-slice-label');
  const coronalCanvas = document.getElementById('coronal-canvas');
  const axialCanvas = document.getElementById('axial-canvas');

  const Z = appState.volume.shape[2];
  if (axialSlider) {
    axialSlider.min = 0;
    axialSlider.max = Z - 1;
    axialSlider.value = appState.currentAxialSlice;
    axialSlider.oninput = () => {
      appState.currentAxialSlice = parseInt(axialSlider.value);
      if (axialLabel) axialLabel.textContent = `Slice ${appState.currentAxialSlice} / ${Z-1}`;
      renderAxialSlice(axialCanvas, appState.volume, appState.mask, appState.currentAxialSlice, appState.showMask);
    };
    if (axialLabel) axialLabel.textContent = `Slice ${appState.currentAxialSlice} / ${Z-1}`;
  }

  // Toggle mask button
  const toggleBtn = document.getElementById('toggle-mask-btn');
  if (toggleBtn) {
    toggleBtn.onclick = () => {
      appState.showMask = !appState.showMask;
      toggleBtn.textContent = appState.showMask ? 'Hide Overlay' : 'Show Overlay';
      toggleBtn.classList.toggle('btn-active', appState.showMask);
      renderAxialSlice(axialCanvas, appState.volume, appState.mask, appState.currentAxialSlice, appState.showMask);
    };
  }

  // Render canvases after a short delay (DOM must be ready)
  setTimeout(() => {
    if (axialCanvas) {
      renderAxialSlice(axialCanvas, appState.volume, appState.mask, appState.currentAxialSlice, appState.showMask);
    }
    if (coronalCanvas && results.callosalSlice >= 0) {
      renderCoronalSlice(coronalCanvas, appState.volume, appState.mask, results.callosalSlice, results.callosalData);
    }
  }, 100);
}

// ─── Sample Data Loader ──────────────────────────────────────────────────────

async function loadSampleData() {
  showScreen('processing');
  const fnEl = document.getElementById('processing-filename');
  if (fnEl) fnEl.textContent = 'Sample CT — CADS BrainCT-1mm Subject 155';

  const steps = [
    'Loading sample data',
    'Building brain mask',
    'Extracting CSF voxels',
    'Morphological filtering',
    'Isolating ventricles',
    'Computing Evans Index',
    'Computing callosal angle',
    'Computing volume',
    'Generating report'
  ];
  initProgress(steps);

  try {
    advanceProgress(0, 'Fetching sample CT scan (~430 KB)...');

    const resp = await fetch('./sample-data.json');
    if (!resp.ok) throw new Error('Failed to load sample data');
    const sample = await resp.json();

    advanceProgress(0, 'Decompressing sample volume...');
    await delay(30);

    // Decode base64 -> gzip -> int16 -> float32
    const b64 = sample.data_b64_gzip_int16;
    const binaryStr = atob(b64);
    const compressed = new Uint8Array(binaryStr.length);
    for (let i = 0; i < binaryStr.length; i++) compressed[i] = binaryStr.charCodeAt(i);

    const raw = pako.inflate(compressed);
    const int16 = new Int16Array(raw.buffer);
    const float32 = new Float32Array(int16.length);
    for (let i = 0; i < int16.length; i++) float32[i] = int16[i];

    const volume = {
      shape: sample.shape,
      spacing: sample.spacing,
      affine: [
        [sample.spacing[0], 0, 0, 0],
        [0, sample.spacing[1], 0, 0],
        [0, 0, sample.spacing[2], 0],
        [0, 0, 0, 1]
      ],
      data: float32,
      header: {
        ndim: 3,
        datatype: 16,
        bitpix: 32,
        voxOffset: 352,
        sformCode: 0,
        dims: sample.shape,
        pixdim: sample.spacing
      }
    };

    appState.volume = volume;
    appState.fileName = 'sample_ct_155.nii.gz';
    appState.fileSize = compressed.length;

    const results = await runPipeline(volume);
    buildResultsUI(results);
    showScreen('results');

  } catch (err) {
    console.error('Sample data error:', err);
    showError(err.message || 'Failed to load sample data.');
  }
}

// ─── Init ─────────────────────────────────────────────────────────────────────

function resetApp() {
  appState.volume = null;
  appState.mask = null;
  appState.results = null;
  appState.currentAxialSlice = 0;
  appState.showMask = true;
  appState.medsam2Box = null;
  appState.segmentationMethod = 'threshold';
  showScreen('upload');
}

document.addEventListener('DOMContentLoaded', () => {
  // File input
  const fileInput = document.getElementById('file-input');
  const dropZone  = document.getElementById('drop-zone');

  if (fileInput) {
    fileInput.addEventListener('change', e => {
      if (e.target.files && e.target.files.length > 0) {
        handleFiles(e.target.files);
      }
    });
  }

  // Drag and drop
  if (dropZone) {
    dropZone.addEventListener('click', () => fileInput && fileInput.click());

    dropZone.addEventListener('dragover', e => {
      e.preventDefault();
      dropZone.classList.add('drag-over');
    });

    dropZone.addEventListener('dragleave', () => {
      dropZone.classList.remove('drag-over');
    });

    dropZone.addEventListener('drop', e => {
      e.preventDefault();
      dropZone.classList.remove('drag-over');
      if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
        handleFiles(e.dataTransfer.files);
      }
    });
  }

  // New scan button
  const newScanBtn = document.getElementById('new-scan-btn');
  if (newScanBtn) newScanBtn.addEventListener('click', resetApp);

  // Sample data button
  const sampleBtn = document.getElementById('sample-btn');
  if (sampleBtn) sampleBtn.addEventListener('click', loadSampleData);

  // MedSAM2 server config
  const serverUrlInput  = document.getElementById('server-url-input');
  const checkServerBtn  = document.getElementById('check-server-btn');

  // Set server URL from the input field on init (now defaults to HF Space)
  if (serverUrlInput && serverUrlInput.value.trim()) {
    MedSAMClient.setServerUrl(serverUrlInput.value.trim());
  }

  if (serverUrlInput) {
    serverUrlInput.addEventListener('change', () => {
      const url = serverUrlInput.value.trim();
      if (url) MedSAMClient.setServerUrl(url);
      updateServerStatusUI(null); // reset status when URL changes
    });
  }

  if (checkServerBtn) {
    checkServerBtn.addEventListener('click', async () => {
      const url = serverUrlInput ? serverUrlInput.value.trim() : '';
      if (url) MedSAMClient.setServerUrl(url);

      checkServerBtn.disabled = true;
      checkServerBtn.textContent = 'Checking…';
      updateServerStatusUI('checking');

      const result = await MedSAMClient.checkHealth();
      appState.medsam2Available = result.available;
      updateServerStatusUI(result.available ? 'connected' : 'disconnected');

      checkServerBtn.disabled = false;
      checkServerBtn.textContent = 'Check Connection';
    });
  }

  // Show upload screen
  showScreen('upload');
});

/**
 * Update the server status indicator in the upload screen.
 * state: 'connected' | 'disconnected' | 'checking' | null (unknown)
 */
function updateServerStatusUI(state) {
  const el = document.getElementById('server-status');
  if (!el) return;
  el.className = 'server-status'; // reset classes
  if (state === 'connected') {
    el.classList.add('server-status-connected');
    el.textContent = 'Connected';
  } else if (state === 'disconnected') {
    el.classList.add('server-status-disconnected');
    el.textContent = 'Not reachable';
  } else if (state === 'checking') {
    el.classList.add('server-status-unknown');
    el.textContent = 'Checking…';
  } else {
    el.classList.add('server-status-unknown');
    el.textContent = 'Not connected';
  }
}
