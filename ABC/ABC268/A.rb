require 'set'
nums = gets.split.map(&:to_i)
s = Set[*nums]
puts s.size()