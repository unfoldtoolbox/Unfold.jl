using LiveServer

dr = filter(isdir,readdir(joinpath("src","generated"),join=true))
push!(dr,"./build")
servedocs(skip_dirs=dr,literate_dir=joinpath("literate"),foldername=".",host="0.0.0.0")


