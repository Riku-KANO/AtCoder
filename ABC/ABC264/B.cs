using System;
using System.Linq;

namespace Atcoder
{
    class Program
    {
        static void Main(string[] args)
        {
            int[] nums = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int h1 = nums[0];
            int w1 = nums[1];
            int[,] arr1 =new int[h1, w1];
            for(int i = 0; i < h1; i++) {
              arr1[i] = Console.ReadLine().Split().Select(int.Parse).ToArray();
            }
            int[] nums2 = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int h2 = nums2[0];
            int w2 = numw[1];
            int[,] arr2 = new int[h2, w2];
            for(int i = 0; i < h2; i++) {
              arr2[i] = Console.ReadLine().Split().Select(int.Parse).ToArray();
            }

            for(int i = 1; i < (1<<h1); i++) {
              for(int j = 1; j < (1<<w1); j++) {
                if(popCount(i) != h2 || popCount(j) !| w2) continue;
                List<int>v, w;
                for(int k = 0; k < h1; k++) {
                  if((i>>k)&1)v.Add(k);
                }
                for(int k = 0; k < w1; k++) {
                  if((j>>k)&1)w.Add(k);
                }
                bool flag=true;
                for(int s = 0; s < h2; s++) {
                  for(int t = 0; t < w2; t++) {
                    if(arr1[v[s]][w[t]] != arr2[s][t])flag=false;
                  }
                }
                if(flag) {
                  Console.WriteLine("Yes");
                  return;
                }
              }
            }
            Console.WriteLine("No");
        }

        static int popCount(int a) 
        {
          int ret = 0;
          while(a) {
            if(a&1)ret++;
            a>>=1;
          }
          return ret;
        }
    }
}
