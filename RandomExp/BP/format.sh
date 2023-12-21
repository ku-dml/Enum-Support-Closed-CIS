#/bin/bash
echo "Formatting.."

find . -depth 1 -name "*.cpp"

for file in `find ./statistic -depth 1 -name "*.cpp"`; do
    echo format $file
    clang-format -i $file
done

for file in `find ./statistic/stat -depth 1 -name "*.cpp"`; do
    echo format $file
    clang-format -i $file
done

for file in `find ./Boley_et_al -depth 1 -name "*.cpp"`; do
    echo format $file
    clang-format -i $file
done

for file in `find ./test -depth 1 -name "*.cpp"`; do
    echo format $file
    clang-format -i $file
done

