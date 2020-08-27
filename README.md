# nf-training-public

Nextflow Training material publically availible. 

## Get started page 

https://www.seqera.io/training/


* Add read permision for the training bucket ie. `s3://seqeralabs.com/public` by modifying 
the IAM policy `nf-training-ro`. E.g.

```
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": "s3:ListBucket",
            "Resource": "arn:aws:s3:::seqeralabs.com"
        },
        {
            "Sid": "VisualEditor1",
            "Effect": "Allow",
            "Action": [
                "s3:GetObjectAcl",
                "s3:GetObject",
                "s3:GetObjectTagging"
            ],
            "Resource": "arn:aws:s3:::seqeralabs.com/public/*"
        }
    ]
}
```


## Course scripts & data 

All data and scripts should be hosted at this bucket: `s3://seqeralabs.com/public/nf-training`. 

1. Sync the local data and scripts with the bucket:

```
aws s3 sync --acl public-read --delete --include "nf-training/{data,0?_*}" nf-training s3://seqeralabs.com/public/nf-training
```

2. Sync training material: 

```
aws s3 sync s3://seqeralabs.com/public/nf-training .
```

## Publish to seqera.io

1. Ensure the latest version of the Seqera website is located at: `../seqera-website`.

2. Go into `../seqera-website`.

3. Run:

```
make publish invalidate
```

Note the aws command may contain a profile. 


## Misc

Remove unnecessary containers:

```
docker rmi $(docker images -q)
```

Resize the instance volume:

```
bash misc/resize.sh
```

