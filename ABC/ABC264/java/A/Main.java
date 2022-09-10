import java.util.*;

public class Main {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    int[]nums = new int[2];
    Arrays.setAll(nums, i -> sc.nextInt());
    int l = nums[0]-1;
    int r = nums[1];
    System.out.println("atcoder".substring(l, r));
  }
}