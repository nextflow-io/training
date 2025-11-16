// params.models = "https://dl.fbaipublicfiles.com/MMPT/metaclip/b32_400m.pt"
params.models = "data/models/b32_400m.pt"
params.width = 400

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    models = channel.fromPath("data/models/b32_400m.pt")
    | map { model -> [[modelId: model.baseName], model] }

    images
    | combine(models)
    | map { imgMeta, img, modelMeta, model -> [imgMeta + modelMeta, model, img]}
    | Classify
    | view

    return

    Resize(images, params.width)
    | view

}

process Resize {
    input:
        tuple val(meta), path(img)
        val(width)
    output: tuple val(meta), path("resized.png")
    script: "convert $img -resize ${width}x resized.png"
}

process Classify {
    memory '8G'
    input: tuple val(meta), path(model), path(pic)
    output: tuple val(meta), path('summary.json')
    script: "classify.py --model-path $model --image-dir . --json > summary.json"
}

process Measure {
    input: path(img)
    output: tuple val(img.baseName), path('dimensions.txt')
    script:
    """
    identify -format "%wx%h" "$img" > dimensions.txt
    """
}

