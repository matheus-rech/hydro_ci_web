/**
 * MedSAM2 API Client
 * HydroMorph — hydrocephalus morphometrics pipeline
 *
 * Communicates with a local or remote MedSAM2 inference server.
 * All segmentation is optional — the app falls back to the threshold
 * pipeline when the server is unavailable.
 *
 * Requires: pako (for gzip compression/decompression)
 */

const MedSAMClient = (() => {

  // Default server URL — user can override via setServerUrl()
  let serverUrl = 'http://localhost:5000';

  /** Update the server base URL (strips trailing slash). */
  function setServerUrl(url) {
    serverUrl = (url || '').trim().replace(/\/$/, '');
  }

  /** Return the current server URL. */
  function getServerUrl() {
    return serverUrl;
  }

  /**
   * Check whether the MedSAM2 server is reachable.
   * Returns { available: boolean, ...serverInfo }
   * Timeout: 5 seconds.
   */
  async function checkHealth() {
    try {
      const resp = await fetch(`${serverUrl}/api/health`, {
        signal: AbortSignal.timeout(5000),
      });
      if (!resp.ok) return { available: false };
      const info = await resp.json().catch(() => ({}));
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

    // ── 1. Window to brain CT soft-tissue range and normalise to 0–255 ────────
    // Default window: C=40 W=80  →  [0, 80] HU
    const WL = 40;   // window level (centre)
    const WW = 80;   // window width
    const loHU = WL - WW / 2;  //  0
    const hiHU = WL + WW / 2;  // 80

    const windowed = new Uint8Array(total);
    for (let i = 0; i < total; i++) {
      const hu = volumeData[i];
      const clamped = Math.max(loHU, Math.min(hiHU, hu));
      windowed[i] = Math.round(((clamped - loHU) / WW) * 255);
    }

    // ── 2. Gzip-compress with pako ────────────────────────────────────────────
    const compressed = pako.gzip(windowed);

    // ── 3. Build FormData ─────────────────────────────────────────────────────
    const formData = new FormData();
    formData.append('volume', new Blob([compressed], { type: 'application/octet-stream' }), 'volume.raw.gz');
    formData.append('shape',   JSON.stringify(shape));
    formData.append('spacing', JSON.stringify(spacing));
    formData.append('box',     JSON.stringify(box));

    // ── 4. POST to server ─────────────────────────────────────────────────────
    const resp = await fetch(`${serverUrl}/api/segment`, {
      method: 'POST',
      body: formData,
    });

    if (!resp.ok) {
      let errMsg = `Server returned ${resp.status}`;
      try {
        const errJson = await resp.json();
        if (errJson && errJson.error) errMsg = errJson.error;
      } catch { /* ignore parse error */ }
      throw new Error(errMsg);
    }

    const result = await resp.json();

    // ── 5. Decode returned mask (base64 → gzip → Uint8Array) ─────────────────
    const maskB64 = result.mask_b64_gzip;
    if (!maskB64) {
      throw new Error('Server response did not include mask_b64_gzip field.');
    }

    const binaryStr     = atob(maskB64);
    const maskCompressed = new Uint8Array(binaryStr.length);
    for (let i = 0; i < binaryStr.length; i++) {
      maskCompressed[i] = binaryStr.charCodeAt(i);
    }

    const mask = pako.inflate(maskCompressed);
    return new Uint8Array(mask);
  }

  return { setServerUrl, getServerUrl, checkHealth, segment };

})();
