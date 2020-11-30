################################################################################  
### Variables 

variable instance_name {
  description = "Name the instance"
  default     = "SANGER_USERNAME-projman"
}

variable instance_image_name {
  description = "Image to boot this instance from"
  default     = "bionic-WTSI-docker_55329_884ace11"
}

variable instance_flavor_name {
  description = "Flavor for this instance"
  default     = "m1.xlarge"
}

variable ssh_user {
  description = "SSH user name"
  default     = "ubuntu"
}

variable private_key {
  description = "Path to private SSH key"
  default     = "keys/user_key"
}

variable public_key {
  description = "Path to private SSH key"
  default     = "keys/user_key.pub"
}

variable openstack_keypair_name {
  description = "Name of the keypair for this machine (must be present in openstack already)"
  default     = "SANGER_USERNAME-projman"
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
resource "openstack_compute_keypair_v2" "user_porjman_keypair" {
  name       = var.openstack_keypair_name
  public_key = file(var.public_key)
}

################################################################################  
### Create instance 
resource "openstack_compute_instance_v2" "user_instance" {
  depends_on = [ openstack_compute_keypair_v2.user_porjman_keypair ]
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
resource "openstack_compute_floatingip_v2" "user_floating_ip" {
  pool  = var.floating_ip_pool
}

################################################################################  
### Asociate Floting IP 
resource "openstack_compute_floatingip_associate_v2" "associate_floating_ip" {
  depends_on = [ openstack_compute_floatingip_v2.user_floating_ip, openstack_compute_instance_v2.user_instance ]
  floating_ip = openstack_compute_floatingip_v2.user_floating_ip.address
  instance_id = openstack_compute_instance_v2.user_instance.id
}

################################################################################  
### Provision the instance with nginx and required  files
resource null_resource "provision_instance" {
  depends_on = [ openstack_compute_floatingip_associate_v2.associate_floating_ip ]

  connection {
    user         = var.ssh_user
    port         = 22
    private_key  = file(var.private_key)
    host         = openstack_compute_floatingip_v2.user_floating_ip.address
  }

  provisioner "file" {
    source     = "${path.root}/provision/projman"
    destination = "/tmp/projman"
  }

  provisioner "file" {
    source     = "${path.root}/provision/provision.sh"
    destination = "/tmp/provision.sh"
  }

  provisioner "remote-exec" {
    inline = [
      "sudo mv /tmp/projman /projman",
      "sudo chown -R $USER:$USER /projman",
      "sudo chmod +x /tmp/provision.sh",
      "sudo /tmp/provision.sh SANGER_USERNAME ${openstack_compute_floatingip_v2.user_floating_ip.address}",
    ]
  }
}

output "user_server_ip" {
  value = openstack_compute_floatingip_v2.user_floating_ip.address
}
output "user_server_name" {
  value = "SANGER_USERNAME"
}
