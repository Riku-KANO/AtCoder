// package ABC.ABC270.java.B;
import java.util.*;

public class Main {
  public static void main(String [] args) {
    Scanner sc = new Scanner(System.in);
    Integer x, y, z;
    x = sc.nextInt();
    y = sc.nextInt();
    z = sc.nextInt();
    if(y < 0) {
      x=-x;
      y=-y;
      z=-z;
    }
    if(x < y) {
      System.out.println(Math.abs(x));
    } else {
      if(z > y) {
        System.out.println(-1);
      } else {
        System.out.println(Math.abs(z) + Math.abs(x-z));
      }
    }
  }  
}
