// Python (Pyodide) in a worker, with gpxiltools, and three ways to give it a file:
//   blob  - a File or Blob (as the browser keeps them) mounted with Emscripten's WORKERFS: reads go to the disk,
//           a slice at a time, and nothing is copied into Python's memory;
//   web   - a file on a web server, read with HTTP range requests (Emscripten's createLazyFile, 1 MB chunks);
//   copy  - the whole file copied into Python's in-memory file system (what most in-browser Pythons do).
importScripts("https://cdn.jsdelivr.net/pyodide/v0.29.5/full/pyodide.js");

let py = null;
const wasmBytes = () => (py ? py._module.HEAP8.buffer.byteLength : 0);

// count what is read from the disk (WORKERFS reads through FileReaderSync) and from the web (XHR ranges)
const counts = { diskReads: 0, diskBytes: 0, webRequests: 0, webBytes: 0 };
const readSlice = FileReaderSync.prototype.readAsArrayBuffer;
FileReaderSync.prototype.readAsArrayBuffer = function (blob) {
  counts.diskReads += 1;
  counts.diskBytes += blob.size;
  return readSlice.call(this, blob);
};
const xhrSend = XMLHttpRequest.prototype.send;
XMLHttpRequest.prototype.send = function (...a) {
  this.addEventListener("load", () => {
    counts.webRequests += 1;
    const r = this.response;
    counts.webBytes += typeof r === "string" ? r.length : r ? r.byteLength || 0 : 0;
  });
  return xhrSend.apply(this, a);
};

const handlers = {
  async boot({ wheels }) {
    const t = performance.now();
    py = await loadPyodide();
    await py.loadPackage(["numpy", "pillow", "jsonschema", "micropip"]);
    const micropip = py.pyimport("micropip");
    // deps=False: gpxiltools asks for imagecodecs, which Pyodide does not have; tifffile reads and writes
    // deflate (the gpxil default) with Python's own zlib
    await micropip.install.callKwargs(wheels.map((w) => new URL(w, self.location).href), { deps: false });
    py.runPython("import gpxiltools, tifffile, numpy");
    return { ms: performance.now() - t, wasm: wasmBytes(), versions: py.runPython("f'gpxiltools {gpxiltools.__version__}, tifffile {tifffile.__version__}, numpy {numpy.__version__}'") };
  },
  mountBlob({ dir, name, blob }) {
    py.FS.mkdirTree(dir);
    py.FS.mount(py.FS.filesystems.WORKERFS, { blobs: [{ name, data: blob }] }, dir);
    return { wasm: wasmBytes() };
  },
  mountUrl({ dir, name, url }) {
    py.FS.mkdirTree(dir);
    py.FS.createLazyFile(dir, name, url, true, false);
    return { wasm: wasmBytes() };
  },
  async copyBlob({ dir, name, blob }) {
    py.FS.mkdirTree(dir);
    py.FS.writeFile(`${dir}/${name}`, new Uint8Array(await blob.arrayBuffer()));
    return { wasm: wasmBytes() };
  },
  async run({ code }) {
    for (const k of Object.keys(counts)) counts[k] = 0;
    const t = performance.now();
    let out = await py.runPythonAsync(code);
    if (out && out.toJs) out = out.toJs({ dict_converter: Object.fromEntries });
    return { ms: performance.now() - t, wasm: wasmBytes(), out, ...counts };
  },
};

self.onmessage = async (e) => {
  const { id, cmd, args } = e.data;
  try {
    postMessage({ id, ok: true, result: await handlers[cmd](args || {}) });
  } catch (err) {
    postMessage({ id, ok: false, error: String(err && err.message ? err.message : err).slice(-2000) });
  }
};
