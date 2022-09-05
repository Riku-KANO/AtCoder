using System;
using System.Linq;

namespace Atcoder
{
    class Program
    {
        static void Main(string[] args)
        {
            var s = "atcoder";
            int[] num = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int n = num[1] - num[0] + 1;
            Console.WriteLine(s.Substring(num[0]-1, n));
        }
    }
}
