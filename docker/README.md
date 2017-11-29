# Everything needed for solveme

Docker Image providing solveme, cobrame, ecolime

## Run the image (to build the image see **Building the Image**)
```
docker run --rm -i -v $(pwd):/workdir -t sbrg/solveme:3.6 bash
```
From here you can run python and follow the instructions for importing and running solveme, on
https://github.com/SBRG/solveme or
https://github.com/SBRG/solvemepy


## Building the Image
1. Copy qminos1114b.tar.gz into this folder (i.e., [solveme_root]/docker/
1. Build the image
```
$ docker build -t sbrg/solveme:3.6 .
```

## Transfer the image to another computer
1. [optional] Save the image
```
docker save -o <image_file> sbrg/solveme:3.6
```
2. Copy (rsync, scp) the <image_file> above to your computer. You might want to tar/gzip it first.
3. Load the image
```
docker load -i <image_file>
```
