#!/bin/zsh

set -u

ATLAS_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_MARKER='name="haa-build" content="v3-internal"'
PORT=""

for candidate in 8793 8794 8795; do
  if curl -fsS --max-time 2 "http://127.0.0.1:${candidate}/" 2>/dev/null |
    grep -Fq "${BUILD_MARKER}"; then
    ATLAS_URL="http://127.0.0.1:${candidate}/?launch=local#/atlas"
    open -a "Google Chrome" "${ATLAS_URL}"
    echo "The current Human Aging Atlas is already running."
    echo "You can close this window."
    exit 0
  fi

  if ! lsof -nP -iTCP:"${candidate}" -sTCP:LISTEN >/dev/null 2>&1; then
    PORT="${candidate}"
    break
  fi
done

if [[ -z "${PORT}" ]]; then
  echo "The Atlas could not find an available local port."
  echo "Close older Atlas Terminal windows and try again."
  exit 1
fi

cd "${ATLAS_DIR}" || exit 1

python3 -m http.server "${PORT}" --bind 127.0.0.1 &
SERVER_PID=$!
ATLAS_URL="http://127.0.0.1:${PORT}/?launch=local#/atlas"

stop_server() {
  kill "${SERVER_PID}" >/dev/null 2>&1 || true
}

trap stop_server EXIT INT TERM

for attempt in {1..20}; do
  if curl -fsS --max-time 1 "http://127.0.0.1:${PORT}/" 2>/dev/null |
    grep -Fq "${BUILD_MARKER}"; then
    open -a "Google Chrome" "${ATLAS_URL}"
    echo ""
    echo "Human Aging Atlas is running:"
    echo "${ATLAS_URL}"
    echo ""
    echo "Keep this Terminal window open while using the Atlas."
    echo "Press Control-C or close this window to stop the local server."
    wait "${SERVER_PID}"
    exit 0
  fi
  sleep 0.25
done

echo "The Atlas could not start on port ${PORT}."
echo "Close this window and double-click the launcher again."
exit 1
