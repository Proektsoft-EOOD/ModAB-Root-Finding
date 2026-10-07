#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
dotnet build tests/ModAB.Switching.Verification -c Release -m:1 -nr:false
dotnet tests/ModAB.Switching.Verification/bin/Release/net10.0/ModAB.Switching.Verification.dll --output results/verification.json
