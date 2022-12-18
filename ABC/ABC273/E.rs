use proconio::{*};
use std::{collections::HashMap};

struct Node {
  val: i32,
  par: usize,
  vec: Vec<i32>,
} 

#[fastout]
fn solve1() {
  input! {
    q: usize
  }
  let mut nodes: Vec<Node> = Vec::new();
  let mut notebook: HashMap<usize, usize> = HashMap::new();
  nodes.push(Node{val: -1, par: 0, vec:Vec::new()});
  let mut node_id: usize = 0;
  for _ in 0..q {
    input! {
      query: String
    }
    if query == "ADD" {
      input! {
        val: i32
      } 
      let new_node = Node{val: val, par: node_id, vec:Vec::new()};
      let new_node_id = nodes.len();
      nodes.push(new_node);
      nodes[node_id].vec.push((new_node_id - 1) as i32);
      node_id = new_node_id;
    }
    else if query == "DELETE" {
      if node_id != 0 {
        node_id = nodes[node_id].par;
      }
    }
    else if query == "SAVE" {
      input! {
        page: usize
      }
      notebook.entry(page).
        and_modify(|e| *e = node_id)
        .or_insert(node_id);

    }
    else if query == "LOAD" {
      input! {
        page: usize
      }
      notebook.entry(page).or_insert(0);
      let new_id: usize = notebook[&page];
      node_id = new_id;
    }
    println!("{}", nodes[node_id as usize].val); 
  }
}

#[fastout]
fn solve2() -> Result<(),()> {
  let mut buf = String::new();
  std::io::stdin().read_line(&mut buf).unwrap();
  let q = buf.trim().parse().unwrap();

  let mut nodes: Vec<Node> = Vec::new();
  let mut notebook: HashMap<usize, usize> = HashMap::new();
  nodes.push(Node{val: -1, par: 0, vec:Vec::new()});
  let mut node_id: usize = 0;

  // query input
  let queries: Vec<Vec<String>> = (0..q).map(|_| {
    let mut buf = String::new();
    std::io::stdin().read_line(&mut buf).unwrap();
    let line: Vec<String> = buf.trim().split_whitespace().map(|s|s.to_string()).collect();
    line
  }).collect();


  // process
  let mut result: Vec<i32> = Vec::new();
  for v in queries.iter() {
    match v[0].as_str() {
      "ADD" => {
        let val: i32 = v[1].parse().unwrap();
        let new_node = Node{val: val, par: node_id, vec:Vec::new()};
        let new_node_id = nodes.len();
        nodes.push(new_node);
        nodes[node_id].vec.push((new_node_id - 1) as i32);
        node_id = new_node_id;
      }, 
      "DELETE" => {
        if node_id != 0 {
          node_id = nodes[node_id].par;
        }
      },
      "SAVE" => {
        let page: usize = v[1].parse().unwrap();
        notebook.entry(page)
          .and_modify(|e| *e = node_id)
          .or_insert(node_id);
      },
      "LOAD" => {
        let page: usize = v[1].parse().unwrap();
        notebook.entry(page).or_insert(0);
        let new_id: usize = notebook[&page];
        node_id = new_id;
      }
      _ =>{ unreachable!() }
    }
    result.push(nodes[node_id].val);
  }
  println!("{}", result.iter().map(|s|s.to_string()).collect::<Vec<String>>().join(" "));
  Ok(())
}

fn main() {
  solve2();
}