using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;

namespace AtCoder
{
    class Program
    {
        static void Main(string[] args)
        {
            int[] nums = Console.ReadLine().Split().Select(int.Parse).ToArray();
            HashSet<int> hs = new HashSet<int>(nums);
            Console.Out.WriteLine(hs.Count());
        }
    }
}
