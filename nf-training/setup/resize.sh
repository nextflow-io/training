#!/bin/bash
#
# Resize Cloud9 boot disk
# https://docs.aws.amazon.com/cloud9/latest/user-guide/move-environment.html
#

# Specify the desired volume size in GiB as a command-line argument. If not specified, default to 20 GiB.
SIZE=${1:-20}

# Install the jq command-line JSON processor.
sudo yum -y install jq

# Get the ID of the envrionment host Amazon EC2 instance.
INSTANCEID=$(curl http://169.254.169.254/latest/meta-data//instance-id)

# Get the region
X_ZONE=$(curl -fs http://169.254.169.254/latest/meta-data/placement/availability-zone)
export AWS_DEFAULT_REGION="`echo \"$X_ZONE\" | sed 's/[a-z]$//'`"

# Get the ID of the Amazon EBS volume associated with the instance.
VOLUMEID=$(aws ec2 describe-instances --instance-id $INSTANCEID | jq -r .Reservations[0].Instances[0].BlockDeviceMappings[0].Ebs.VolumeId)

# Resize the EBS volume.
aws ec2 modify-volume --volume-id $VOLUMEID --size $SIZE

# Wait for the resize to finish.
while [ "$(aws ec2 describe-volumes-modifications --volume-id $VOLUMEID --filters Name=modification-state,Values="optimizing","completed" | jq '.VolumesModifications | length')" != "1" ]; do
  sleep 1
done

# Rewrite the partition table so that the partition takes up all the space that it can.
sudo growpart /dev/xvda 1

# Expand the size of the file system.
sudo resize2fs /dev/xvda1
