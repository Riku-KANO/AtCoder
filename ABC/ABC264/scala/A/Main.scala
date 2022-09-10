import java.util.Scanner

object Main {
  def main(args: Array[String]): Unit = {
    val sc = new Scanner(System.in)
    val l = sc.nextInt - 1
    val r = sc.nextInt
    println("atcoder".substring(l, r))
  }
}
