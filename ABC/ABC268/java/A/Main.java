import java.util.*;

class Main {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    Set<Integer> s = new HashSet<Integer>();
    for(Integer i = 0; i < 5; i++) {
      Integer a = sc.nextInt();
      s.add(a);
    }
    System.out.println(s.size());
  }
}