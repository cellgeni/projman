#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "📡 pre-provision for proxy"
if [[ -d "$DIR/proxy" ]]; then
    echo "✋ proxy folder already exists, aborting..."
    exit 1
else
    cp -r "$DIR/.templates/proxy-server/" "$DIR/proxy/"
fi

cd "$DIR/proxy"
./create_proxy_instance.sh

cd "$DIR"
