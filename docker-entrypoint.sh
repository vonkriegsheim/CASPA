#!/usr/bin/env bash
# CASPA container entrypoint. Dispatches on the first argument:
#   gui     -> launch the Streamlit GUI on :8501 (bound to 0.0.0.0 for the host)
#   doctor  -> run the install checker
#   <else>  -> run the pipeline (default), e.g. --workdir /work --cores 8
set -e
case "$1" in
  gui)
    shift
    exec streamlit run /opt/caspa/caspa/gui.py \
         --server.address 0.0.0.0 --server.port 8501 --server.headless true "$@"
    ;;
  doctor)
    shift
    exec python /opt/caspa/caspa/doctor.py "$@"
    ;;
  *)
    exec python /opt/caspa/caspa/run.py "$@"
    ;;
esac
