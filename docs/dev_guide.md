# Development guidelines:

For the docs, we use [MkDocs](https://www.mkdocs.org/) because of its flexibility:  
- [mkdocs-material](https://squidfunk.github.io/mkdocs-material/getting-started/): look and feel  
- [mkdocs-bibtex](https://github.com/shyamd/mkdocs-bibtex): literature reference  
- [MkPDFs](https://comwes.github.io/mkpdfs-mkdocs-plugin/getting-started.html): PDF version  


## Setting up mkdocs 

```bash
# osx specific settings
conda install pango cairo

pip install mkdocs
pip install mkdocs-material
pip install mkdocs-bibtex
pip install -e git+https://github.com/jwaschkau/mkpdfs-mkdocs-plugin.git#egg=mkpdfs-mkdocs-plugin

# osx specific settings
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

## Updating docker containers

The dockerfile for container are located at `docker/` directory. Some of the enviroments use conda recipes, which are located at the `envs/` directory. After updating one recipes, to build and upload the container one should use:

```bash
cd Baltica/
docker build -f <dockerfile> --tag <tag> .
docker push <tag>
```

For example,
- of docker file: `docker/baltica/1.0/Dockerfile`
- tag: `tbrittoborges/baltica:latest`

The containers are hosted at the [docker hub](https://hub.docker.com/repository/docker/tbrittoborges/) and this location can be changed at the 
*container* directive at the snakefiles.  

## Contributing to the documentation

### Modify any of the doc files
```bash
vi docs/setup.md 
```

### Test the changes locally
```bash
mkdocs serve
```

If everything looks fine you can submit a patch or pull-request.

### Deploy changes
This requires permissions from the GitHub organization.
```bash
mkdocs gh-deploy
```


## Testing Baltica 
Baltica's continuous integration testing suite is under development.

