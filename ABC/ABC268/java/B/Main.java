import java.util.*;

class Main {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    String s = sc.nextLine();
    String t = sc.nextLine();
    if(s.length() > t.length()) {
      System.out.println("No");
      return;
    }
    if(t.substring(0, s.length()).matches(s)) {
      System.out.println("Yes");
    } else {
      System.out.println("No");
    }
  }
}
