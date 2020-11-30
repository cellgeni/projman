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

1. Clone this repository 
2. `cd ~/projman`

### Proxy server

1. `./create_proxy.sh`

2. Wait for it to finish.

3. Check the output for the proxy IP.
```
proxy_server_ip = <proxy_ip>
```

4. Test connection with `cd proxy` and then `./connect.sh` or `ssh ubuntu@<proxy_ip> -i keys/projman_key`

### Createa a user server

Creatae new OpenStack instance and provision with necessary scripts.

From the checkout folder run:
`./create_user.sh <user>`

After it finishes **change to the newly created dir**: `cd users/<user>`

Verify output with `terraform output`
```
user_server_ip = <server_ip>
user_server_name = <user>
```

Check connection with the new machine `./connect.sh` or `ssh ubuntu@<server_ip> -i keys/user_key`

Pack keys and ssh instructions in a tar to share: `./share.sh`

#### Customizing user server

Edit `.templates/user-server/main.tf` file. The first section of the file has the varaibles that can be edited.

- `instance_name`: name for the instance, default is "<username>-projman"
- `instance_image_name`: image name present in openstack, default it base "bionic-WTSI-docker_55329_884ace11
- `instance_flavor_name`: flavor for the user instance", default is "m1.xlarge"
- `network_name`: name of the network the user instnace will connect to, default is "cloudforms_network"

### Add user server to proxy entries

`create_user.sh` script automatically adds the new user to the proxy server, however, if you want to add it manuall you can follow this steps: 
1. `cd proxy`
2. `./connect.sh`
3. `./projman/add_user.sh <username> <server_ip>`

### Teardown

To tear down the proxy or a any created user server, change directory to that folder and then run `terraform destry` and input `yes` when prompted.

## porjMan

Run `projMan` from the user server and follow the instructions.
