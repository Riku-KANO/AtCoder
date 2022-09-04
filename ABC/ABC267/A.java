import java.util.*;

public class Main {
   public static void main(String args[]) {
     Scanner scanner = new Scanner(System.in);
     String s = scanner.nextLine();
     List list = new ArrayList<>(Arrays.asList("Friday","Thursday","Wednesday","Tuesday","Monday"));
     int ans;
     ans = list.indexOf(s) + 1;
     System.out.println(ans);
   }
}