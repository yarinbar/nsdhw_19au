# nsd hw3 - ianre657

## Docker image

`docker pull ianre657/nsd-env`


## Running project

1. `cd hw3/ianre657`

2. (in fish shell)

    ```fish
    docker run --rm -it -v (pwd):/nsd ianre657/nsd-env
    ```

3. running test
    + `make test`
    + `../validate.py`
