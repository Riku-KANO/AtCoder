l1, r1, l2, r2 = gets.split.map(&:to_i)
puts [((l1..r1).to_a & (l2..r2).to_a).size - 1, 0].max