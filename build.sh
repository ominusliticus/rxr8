set -x

cwd=`pwd`
src_dir=${cwd}/ReactionNetwork
build_dir=${cwd}/build

if [ -d $build_dir ]; then
    rm -rf $build_dir
fi

mkdir -p $build_dir
cd $build_dir

src_files=`find ${src_dir} -name *.cpp | grep -v main`
echo $src_files

flags="-Wall -Wpedantic -Wextra -std=c++20 -g"
compiler=clang++

for src in $src_files;
do
    $compiler $flags -c $src
done

build_artifacts=`ls -1`
$compiler $flags -o main $build_artifacts ${src_dir}/main.cpp 