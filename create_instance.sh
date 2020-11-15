#!/usr/bin/env bash

CURRENT_DIR="$(pwd)"
USERNAME=$1

if [[ -z "$USERNAME" ]]; then
    echo "👤 Need <username>. Usage: ./create_instance.sh <username>"
    exit 1
fi

echo "🗂️  pre-provision for '$USERNAME'"
if [[ -d "$USERNAME" ]]; then
    echo "❌ folder '$USERNAME' already exists, aborting..."
    exit 1
else
    cp -r user-server/ "$USERNAME/"
fi

cd "$CURRENT_DIR/$USERNAME/"
echo "🚀 starting instance creation process for '$USERNAME'"
./create_instance.sh "$USERNAME"

user_ip=$(terraform output | grep server_ip | awk -F' = ' '{print $2}')
user_name=$USERNAME

cd "$CURRENT_DIR/proxy-server"
./add_user.sh "$user_name" "$user_ip" 

cd "$CURRENT_DIR"
