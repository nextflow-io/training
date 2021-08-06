#!/bin/bash
# 
# Copyright (c) 2019-2021, Seqera Labs.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. 
# 
# This Source Code Form is "Incompatible With Secondary Licenses", as
# defined by the Mozilla Public License, v. 2.0.
# 
set -e
set -u

X_TYPE=t2.medium
X_REGION=eu-central-1
X_EVAN=arn:aws:iam::195996028523:user/evan
X_PAOLO=arn:aws:iam::195996028523:user/yo
X_PWD=Secret123!

if [ $# -eq 0 ]; then
    echo "No arguments supplied"
    exit 1
fi


export AWS_DEFAULT_OUTPUT="text"

#
# create IAM user 
#
create_user() {
    # create user
    X_NAME=$1
    X_USER=$(aws iam create-user --user-name $X_NAME | cut -f 2)
    echo "$X_NAME"

    sleep 5
    # Allow login
    aws iam create-login-profile --user-name $X_NAME --password $X_PWD
    sleep 5

    # add to group
    aws iam add-user-to-group --user-name "$X_NAME" --group-name MyCloud9Group
    sleep 5
}

#
# create Cloud9 env 
# 
function create_env() {
    X_NAME=$1

    # Create Cloud9 env
    X_ENV=$(aws --region $X_REGION cloud9 create-environment-ec2 \
    --name "$X_NAME" \
    --instance-type $X_TYPE \
    --owner-arn arn:aws:iam::195996028523:user/$X_NAME \
    --automatic-stop-time-minutes 60)
  sleep 5 

    # Add user to env
    if [[ $X_NAME != evan ]]; then
    aws --region $X_REGION cloud9 create-environment-membership --environment-id $X_ENV --user-arn $X_EVAN --permission read-write
    sleep 5
    fi
    if [[ $X_NAME != yo ]]; then
    aws --region $X_REGION cloud9 create-environment-membership --environment-id $X_ENV --user-arn $X_PAOLO --permission read-write
    sleep 5
    fi 

    echo "User: $X_NAME; Account: 195996028523; Cloud9_environment: https://${X_REGION}.console.aws.amazon.com/cloud9/ide/$X_ENV"
}

create_all() {
    create_user $1
    create_env $1
}

if [[ $1 == all ]]; then
  for x in $(cat user-names.txt); do 
    create_all $x
  done
elif [[ $1 == yo ]]; then
  X_REGION=eu-west-1
  create_env yo
elif [[ $1 == evan ]]; then
  X_REGION=eu-west-1
  create_env evan
else 
  create_env $1
fi
