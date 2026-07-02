#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
VERBOSE=false
BUILD_DIR="$DIR/build"
project=GROOT

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform a full build
    --include-orion=<path>  Set external ORION path
    --include-finer=<path>  Set external FiNeR path
    --compilers=<name>      Set compilers (intel, gnu)
    --use-openmp            Use OpenMP
    --use-radlib            Enable the SLW model via radlib (lib/third_party/radlib-master)
    --debug                 Build in Debug mode (default: Release)

  compile                   Compile the program using the CMakePresets file

  update                    Download git submodules
    --remote                Use the latest remote commit

EOF
    exit 1
}

log() {
    if [ "$VERBOSE" = true ]; then
        echo -e "\033[1;90m$1\033[0m"
    fi
}

error() {
    echo -e "\033[1;31m[ERROR] $1\033[0m" >&2
}

task() {
    echo -e "\033[1;38;5;186m==> $1\033[0m"
}


# Create default CMakePresets.json after a successful build
function write_presets() {
  FC=$(grep  '^CMAKE_Fortran_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep  '^CMAKE_C_COMPILER:FILEPATH='       "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep '^CMAKE_CXX_COMPILER:FILEPATH='     "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)

  cat <<EOF > CMakePresets.json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 22
  },
  "configurePresets": [
    {
      "name": "default",
      "description": "Default preset",
      "binaryDir": "\${sourceDir}/build",
      "cacheVariables": {
        "ORION_PATH":              "${ORION_PATH}",
        "FINER_PATH":              "${FINER_PATH}",
        "CMAKE_BUILD_TYPE":        "${BUILD_TYPE}",
        "CMAKE_Fortran_COMPILER":  "${FC}",
        "CMAKE_C_COMPILER":        "${CC}",
        "CMAKE_CXX_COMPILER":      "${CXX}",
        "USE_OPENMP":              "${USE_OPENMP}",
        "USE_MPI":                 "${USE_MPI}",
        "USE_RADLIB":              "${USE_RADLIB}"
      }
    }
  ]
}
EOF
}


# ---------------------------------------------------------------------------
# Defaults
COMMAND=""
COMPILERS=""
ORION_PATH=$(pwd)'/lib/ORION/'
FINER_PATH=$(pwd)'/lib/third_party/FiNeR/'
USE_OPENMP="false"
USE_MPI="false"
USE_RADLIB="false"
REMOTE="false"
BUILD_TYPE="RELEASE"

CMD=("build" "compile" "update")
CMD_OPTIONS_build=("--compilers --include-orion --include-finer --use-openmp --use-radlib --debug")
CMD_OPTIONS_update=("--remote")

# Parse global options
while getopts "hv-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                help)    usage ;;
                *)       error "Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        h) usage ;;
        v) VERBOSE=true ;;
        ?) error "Unknown global option '-$OPTARG'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
    error "No command provided!"
    usage
fi

COMMAND="$1"
if [[ ! " ${CMD[@]} " =~ " ${COMMAND} " ]]; then
    error "Unknown command '$COMMAND'"
    usage
fi
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-orion=*)
            [[ "$COMMAND" == "build" ]] || { error "--include-orion is only valid for 'build'"; exit 1; }
            ORION_PATH="${1#*=}"
            ;;
        --include-finer=*)
            [[ "$COMMAND" == "build" ]] || { error "--include-finer is only valid for 'build'"; exit 1; }
            FINER_PATH="${1#*=}"
            ;;
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error "--compilers is only valid for 'build'"; exit 1; }
            if [[ ! "$1" =~ ^--compilers=(intel|gnu)$ ]]; then
                error "Invalid value for --compilers. Valid values are 'intel' or 'gnu'."
                exit 1
            fi
            COMPILERS="${1#*=}"
            ;;
        --use-openmp)
            [[ "$COMMAND" == "build" ]] || { error "--use-openmp is only valid for 'build'"; exit 1; }
            USE_OPENMP="true"
            ;;
        --use-radlib)
            [[ "$COMMAND" == "build" ]] || { error "--use-radlib is only valid for 'build'"; exit 1; }
            USE_RADLIB="true"
            ;;
        --debug)
            [[ "$COMMAND" == "build" ]] || { error "--debug is only valid for 'build'"; exit 1; }
            BUILD_TYPE="DEBUG"
            ;;
        --remote)
            [[ "$COMMAND" == "update" ]] || { error "--remote is only valid for 'update'"; exit 1; }
            REMOTE="true"
            ;;
        *)
            eval "opts=(\"\${CMD_OPTIONS_${COMMAND}[@]}\")"
            error "Unknown option '$1' for command '$COMMAND'. Valid options: ${opts[@]}"
            usage
            ;;
    esac
    shift
done


# ---------------------------------------------------------------------------
# Execute the selected command
case "$COMMAND" in
    build)
        task "Building $project"

        task "Cloning submodules"
        [[ $ORION_PATH == $(pwd)'/lib/ORION/' ]]              && git submodule update --init lib/ORION
        [[ $FINER_PATH == $(pwd)'/lib/third_party/FiNeR/' ]]  && git submodule update --init --recursive lib/third_party/FiNeR
        # radlib is optional (update=none in .gitmodules): only pulled on demand.
        # --checkout is required to override update=none for this explicit fetch.
        [[ "$USE_RADLIB" == "true" ]]                         && git submodule update --init --checkout lib/third_party/radlib

        task "Configuring and building $project"
        if [[ $COMPILERS == "intel" ]]; then
          export FC="ifx"
          export CC="icx"
          export CXX="icpx"
        elif [[ $COMPILERS == "gnu" ]]; then
          export FC="gfortran"
          export CC="gcc"
          export CXX="g++"
        fi

        log "Build dir:   $BUILD_DIR"
        log "Build type:  $BUILD_TYPE"
        log "ORION path:  $ORION_PATH"
        log "FINER path:  $FINER_PATH"
        log "Use OpenMP:  $USE_OPENMP"
        if [[ -z "${FC+x}" || -z "${CXX+x}" || -z "${CC+x}" ]]; then
          log "Compilers not set — CMake will decide."
        else
          log "Compilers: FC=$FC, CXX=$CXX, CC=$CC"
        fi

        rm -rf "$BUILD_DIR"
        cmake -B "$BUILD_DIR" \
              -DORION_PATH="$ORION_PATH" \
              -DFINER_PATH="$FINER_PATH" \
              -DUSE_OPENMP="$USE_OPENMP" \
              -DUSE_MPI="$USE_MPI"       \
              -DUSE_RADLIB="$USE_RADLIB" \
              -DUSE_TECIO=OFF            \
              -DCMAKE_BUILD_TYPE="$BUILD_TYPE" || exit 1
        cmake --build "$BUILD_DIR" || exit 1
        log "[OK] Compilation successful"

        task "Writing CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;

    compile)
        task "Compiling $project using CMakePresets"
        cmake --preset default    || exit 1
        cmake --build "$BUILD_DIR" || exit 1
        log "[OK] Compilation successful"
        ;;

    update)
        task "Updating git submodules"
        if [[ "$REMOTE" == "true" ]]; then
          log "Updating submodules to latest remote commit"
          git submodule update --init --remote
        else
          log "Updating submodules to pinned commit"
          git submodule update --init
        fi
        log "[OK] Submodules updated"
        ;;

    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac
