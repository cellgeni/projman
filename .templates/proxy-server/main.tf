
################################################################################  
### Variables 

variable instance_name {
  description = "Name the instance"
  default     = "projman-proxy"
}

variable instance_image_name {
  description = "Image to boot this instance from"
  default     = "bionic-server"
}

variable instance_flavor_name {
  description = "Flavor for this instance"
  default     = "o2.medium"
}

variable ssh_user {
  description = "SSH user name"
  default     = "ubuntu"
}

variable private_key {
  description = "Path to private SSH key"
  default     = "keys/projman_key"
}

variable public_key {
  description = "Path to private SSH key"
  default     = "keys/projman_key.pub"
}

variable openstack_keypair_name {
  description = "Name of the keypair for this machine (must be present in openstack already)"
  default     = "projman"
}

variable network_name {
  description = "Name of the network to attach this instance to"
  default     = "cloudforms_network"
}
variable floating_ip_pool {
  description = "Name of the floating IP pool"
  default     = "public"
}

################################################################################  
### Add key to openstack
resource "openstack_compute_keypair_v2" "porjman_keypair" {
  name       = var.openstack_keypair_name
  public_key = file(var.public_key)
}

################################################################################  
### Create instance 
resource "openstack_compute_instance_v2" "proxy_instance" {
  depends_on = [ openstack_compute_keypair_v2.porjman_keypair ]
  name            = var.instance_name
  image_name      = var.instance_image_name
  flavor_name     = var.instance_flavor_name
  key_pair        = var.openstack_keypair_name
  security_groups = ["cloudforms_ext_in"]
  network {
    name = var.network_name
  }
}

################################################################################  
### Create Floating IP
resource "openstack_compute_floatingip_v2" "proxy_floating_ip" {
  pool  = var.floating_ip_pool
}

################################################################################  
### Asociate Floting IP 
resource "openstack_compute_floatingip_associate_v2" "associate_floating_ip" {
  depends_on = [ openstack_compute_floatingip_v2.proxy_floating_ip, openstack_compute_instance_v2.proxy_instance ]
  floating_ip = openstack_compute_floatingip_v2.proxy_floating_ip.address
  instance_id = openstack_compute_instance_v2.proxy_instance.id
}

################################################################################  
### Provision the instance with nginx and required  files
resource null_resource "provision_instance" {
  depends_on = [ openstack_compute_floatingip_associate_v2.associate_floating_ip ]

  connection {
    user         = var.ssh_user
    port         = 22
    private_key  = file(var.private_key)
    host         = openstack_compute_floatingip_v2.proxy_floating_ip.address
  }

  provisioner "file" {
    source     = "${path.root}/provision"
    destination = "/home/ubuntu/projman/"
  }

  provisioner "remote-exec" {
    inline = [
      "chmod +x /home/ubuntu/projman/provision.sh",
      "sudo /home/ubuntu/projman/provision.sh",
      "chmod +x /home/ubuntu/projman/add_user.sh"
    ]
  }
}

output "proxy_server_ip" {
  value = openstack_compute_floatingip_v2.proxy_floating_ip.address
}