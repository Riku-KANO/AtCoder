import java.util.Scanner

object Main {
  def main(args: Array[String]) {
    val sc = new Scanner(System.in) 
    val X = sc.nextInt
    val Y = sc.nextInt
    val N = sc.nextInt
    if (X * 3 >= Y) {
      println(Y*(N/3) + X*(N%3))
    } else {
      println(X*N)
    }
  }
}