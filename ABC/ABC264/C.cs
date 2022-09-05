using System;
using System.Linq;
using System.Collections.Generic;

namespace Atcoder
{
    class Program
    {
        static void Main(string[] args)
        {
            int[] nums = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int h1 = nums[0];
            int w1 = nums[1];
            int[,] arr1 = new int[h1, w1];
            for(int i = 0; i < h1; i++) {
              int[] hoge = Console.ReadLine().Split().Select(int.Parse).ToArray();
              for(int j = 0; j < w1; j++) arr1[i,j] = hoge[j];
            }
            int[] nums2 = Console.ReadLine().Split().Select(int.Parse).ToArray();
            int h2 = nums2[0];
            int w2 = nums2[1];
            int[,] arr2 = new int[h2, w2];
            for(int i = 0; i < h2; i++) {
              int[] hoge = Console.ReadLine().Split().Select(int.Parse).ToArray();
              for(int j = 0; j < w2; j++) arr2[i, j] = hoge[j];
            }

            for(int i = 1; i < (1<<h1); i++) {
              for(int j = 1; j < (1<<w1); j++) {
                if(popCount(i) != h2 || popCount(j) != w2) continue;
                List<int>v = new List<int>();
                List<int>w = new List<int>();
                for(int k = 0; k < h1; k++) {
                  if((i >> k) % 2 == 1) v.Add(k);
                }
                for(int k = 0; k < w1; k++) {
                  if((j >> k) % 2 == 1) w.Add(k);
                }
                bool flag=true;
                for(int s = 0; s < h2; s++) {
                  for(int t = 0; t < w2; t++) {
                    if(arr1[v[s], w[t]] != arr2[s, t]) flag=false;
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
          while(a != 0) {
            if(a%2==1)ret++;
            a>>=1;
          }
          return ret;
        }
    }
}
