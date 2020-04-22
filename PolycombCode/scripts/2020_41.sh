cd /polycomb
num=2020_41
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 1.0 4.0 0.1 2020_41
python3 upload.py polycomb $num results/$num
