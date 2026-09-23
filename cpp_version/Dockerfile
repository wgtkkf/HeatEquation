# Official GCC13 docker image
FROM gcc:13

# Install CMake
RUN apt-get update && apt-get install -y cmake

# Create a working directory
WORKDIR /app

# Copy all from local to the working directory
COPY . /app

# Run CMake build in a docker container
RUN mkdir build && cd build && cmake .. && make

# Execute executable, which is defined in the cmake file
CMD ["./build/bmt"]
