#!/usr/bin/env bash

USERNAME=$1

if [[ -z "$USERNAME" ]]; then
    echo "👤 Need <username>. Usage: ./create_instance.sh <username>"
    exit 1
fi


# double check you've got terraform
if [[ -z "$(which terraform)" ]]; then
    echo "🔎 looking for terraform"
    echo "📥 terrform not found, downloading..."
    wget --quiet https://releases.hashicorp.com/terraform/0.13.5/terraform_0.13.5_linux_amd64.zip -O /tmp/terraform.zip
    if [[ -z "$(which unzip)" ]]; then
        sudo apt-get install unzip
    fi
    echo "🏗️ installing terrform"
    sudo unzip /tmp/terraform.zip -d /usr/local/bin/
fi


# generated keys
if [[ ! -f "keys/user_key" ]]; then
    mkdir -p keys
    echo "🔐 generating keys for '$USERNAME'"
    ssh-keygen -b 2048 -t rsa -f keys/user_key -q -N ""
else
    echo "🔑 found existing credentials!"
fi


# source openstack credentials
if [[ -z "${OS_USERNAME}" ]] && [[ -z "${OS_PROJECT_NAME}" ]] && [[ -z "${OS_PASSWORD}" ]] && [[ -z "${OS_AUTH_URL}" ]] && [[ -z "${OS_REGION_NAME}" ]]; then
    echo "🤔 missing openstack credentials in environment"
    echo "🔐 please source openrc.sh file first"
    echo "🌐 go to https://theta.internal.sanger.ac.uk/project/api_access/ and 'Download OpenStack RC File'"
    exit 1
fi


# set username
echo "📝 updating main.tf for '$USERNAME'"
sed -i "s@SANGER_USERNAME@$USERNAME@g" main.tf


# run terraform
echo "🚜 terraform init/validate/plan"
if terraform init 1>/dev/null && terraform validate 1>/dev/null && terraform plan 1>/dev/null; then
    echo "🦄 terraform apply"
    terraform apply -auto-approve
else
    echo "🤔 something went wrong!"
fi
