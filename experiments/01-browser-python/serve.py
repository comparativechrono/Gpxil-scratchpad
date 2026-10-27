"""A static web server with HTTP range requests and cross-origin isolation, for the experiments.

    python serve.py --data /path/to/gpxil/files      # http://localhost:8010/ ; the files at /data/...

Range requests are what lazy files read with; the isolation headers (COOP and COEP) let the page measure the
memory of itself and its workers (performance.measureUserAgentSpecificMemory).
"""

from __future__ import annotations

import argparse
import os
import re
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

_RANGE = re.compile(r"bytes=(\d*)-(\d*)$")


class Handler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map, ".js": "text/javascript", ".tif": "image/tiff", ".whl": "application/zip"}
    data_dir = None
    isolate = True

    def translate_path(self, path):
        if self.data_dir and path.startswith("/data/"):
            rel = path.split("?", 1)[0][len("/data/") :]
            return os.path.join(self.data_dir, *[p for p in rel.split("/") if p not in ("", ".", "..")])
        return super().translate_path(path)

    def end_headers(self):
        if self.isolate:
            self.send_header("Cross-Origin-Opener-Policy", "same-origin")
            self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
            self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def send_head(self):
        rng = self.headers.get("Range")
        path = self.translate_path(self.path)
        if not rng or not os.path.isfile(path):
            return super().send_head()
        m = _RANGE.match(rng.strip())
        size = os.path.getsize(path)
        if not m:
            self.send_error(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            return None
        a, b = m.groups()
        start, end = (max(0, size - int(b)), size - 1) if a == "" else (int(a), min(size - 1, int(b) if b else size - 1))
        if start >= size or start > end:
            self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return None
        fh = open(path, "rb")
        fh.seek(start)
        self.send_response(HTTPStatus.PARTIAL_CONTENT)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.end_headers()
        self._remaining = end - start + 1
        return fh

    def copyfile(self, source, outputfile):
        left = getattr(self, "_remaining", None)
        if left is None:
            return super().copyfile(source, outputfile)
        while left > 0:
            buf = source.read(min(left, 1 << 16))
            if not buf:
                break
            outputfile.write(buf)
            left -= len(buf)
        self._remaining = None

    def log_message(self, *a):
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8010)
    ap.add_argument("--data", help="a folder served at /data/")
    ap.add_argument("--no-isolation", action="store_true")
    a = ap.parse_args()
    Handler.data_dir = a.data
    Handler.isolate = not a.no_isolation
    here = os.path.dirname(os.path.abspath(__file__))
    httpd = ThreadingHTTPServer(("127.0.0.1", a.port), partial(Handler, directory=here))
    print(f"serving {here} on http://localhost:{a.port}/ (data: {a.data})", flush=True)
    httpd.serve_forever()


if __name__ == "__main__":
    main()
