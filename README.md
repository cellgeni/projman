# JupyterHub project launcher

## Overview of the architecture

```ascii
           +-------+
           | proxy |
           +-------+
               +
               |                        +------+
      +--------+---------+        +---> | 🐳 1 |
      |        |         |        |     +------+
      V        V         V        |
  +------+  +-----+  +------+     |     +------+
  | 👤 A |  | 👤 B |  | 👤 C |+----+---> | 🐳 2 |
  +------+  +-----+  +------+     |     +------+
                                  |
                                  |     +------+
                                  +---> | 🐳 3 |
                                        +------+
```

## Provisioning notes

### Proxy server

1. `cd proxy-server`
2. `./create_proxy_instance.sh`

Wait for it to finish.

Verify output with `terraform output`
```
proxy_server_ip = <proxy_ip>
```

 Test connection with `./connect.sh` or `ssh ubuntu@<proxy_ip> -i keys/projman_key`

### Createa a user server

Creatae new OpenStack instance and provision with necessary scripts.

From the checkout folder run:
`./create_intance.sh <user>`

_**Note**: Don't run the script directly from "user-server" folder, that's the template and should not be changed._

After it finishes **change to the newly created dir**: `cd <user>`

Verify output with `terraform output`
```
user_server_ip = <server_ip>
user_server_name = <user>
```

Check connection with the new machine `connect.sh` or
`ssh ubuntu@<server_ip> -i keys/user_key`

Pack keys and ssh instructions in a tar to share: `share.sh`

#### Customizing user server

Edit `main.tf` file. The first section of the file has the varaibles that can be edited.

- `instance_name`: name for the instance, default is "<username>-projman"
- `instance_image_name`: image name present in openstack, default it base "bionic-WTSI-docker_55329_884ace11
- `instance_flavor_name`: flavor for the user instance", default is "m1.xlarge"
- `network_name`: name of the network the user instnace will connect to, default is "cloudforms_network"

### Add user server to proxy entries

When `create_intance.sh` from the root folder is used, the user is added automatically to the proxy server. If you want to add it manuall, you can follow this steps: 
1. `cd proxy-server`
2. `./connect.sh`
3. `./projman/add_user.sh <username> <server_ip>`

### Teardown

To tear down the proxy or a any created user server, change directory to that folder and then run `terraform destry`.

## porjMan

Run `projMan` from the user server and follow the instructions.
