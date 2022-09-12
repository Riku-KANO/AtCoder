import java.util.*;

public class Main {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    int X = sc.nextInt();
    int Y = sc.nextInt();
    int N = sc.nextInt();
    if(X * 3 >= Y) {
      System.out.println((N / 3) * Y + X*(N%3));
    } else {
      System.out.println(X * N);
    }
  }
}