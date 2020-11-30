#!/usr/bin/env bash

CURRENT_DIR="$(pwd)"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
USERNAME=$1

if [[ -z "$USERNAME" ]]; then
    echo "👤 Need <username>. Usage: ./create_instance.sh <username>"
    exit 1
fi

mkdir -p "$DIR/users"

echo "📂 pre-provision for '$USERNAME'"
if [[ -d "$DIR/users/$USERNAME" ]]; then
    echo "✋ folder '$USERNAME' already exists, aborting..."
    exit 1
else
    cp -r "$DIR/.templates/user-server/" "$DIR/users/$USERNAME/"
fi

cd "$DIR/users/$USERNAME/"
echo "🚀 starting instance creation process for '$USERNAME'"
./create_instance.sh "$USERNAME"

user_ip=$(terraform output | grep server_ip | awk -F' = ' '{print $2}')
user_name=$USERNAME

if [[ -d "$DIR/proxy" ]]; then
    cd "$DIR/proxy"
    ./add_user.sh "$user_name" "$user_ip"
else
    echo "✋ proxy folder not found, have you the proxy?"
    cd "$DIR"
    exit 1
fi

cd "$DIR"
