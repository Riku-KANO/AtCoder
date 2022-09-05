using System;
using System.Linq;

namespace Atcoder
{
    class Program
    {
        static void Main(string[] args)
        {
            int[] nums = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int r = nums[0];
            int c = nums[1];
            if(Math.Max(Math.Abs(8-r), Math.Abs(8-c))%2 == 0)
            {
                Console.WriteLine("white");
            } else {
                Console.WriteLine("black");
            }
        }
    }
}
