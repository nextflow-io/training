export AWS_DEFAULT_OUTPUT="text"

confirm() {
    read -r -p "Please confirm you want to destroy ALL Cloud9 environments? [y/N] " response
    case "$response" in
        [yY]|[yY])
        echo ''
        ;;
        *)
        exit 1
        ;;
    esac

}

env_list() {
  aws --region $1 cloud9 list-environments
}

env_destroy() {
  for x in $(env_list $1 | cut -f 2); do 
    aws --region $1 cloud9 delete-environment --environment-id $x; 
  done
}

echo "Will delete the followings Cloud9 envs!"
env_list eu-west-1
env_list eu-central-1
confirm

env_destroy eu-west-1
env_destroy eu-central-1

# echo "Will delete the following user"
# aws iam get-group --group-name MyCloud9Group | grep ^USERS | cut -f 6
# confirm

# for x in $(aws iam get-group --group-name MyCloud9Group | grep ^USERS | cut -f 6); do
#   #aws iam delete-login-profile --user-name $x
#   aws iam detach-group-policy --group-name MyCloud9Group --policy-arn arn:aws:iam::aws:policy/IAMUserChangePassword
#   aws iam delete-user --user-name $x
# done
