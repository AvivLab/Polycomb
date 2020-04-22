cd /polycomb
num=2020_158
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 2.0 1.0 0.2 2020_158
python3 upload.py polycomb $num results/$num
