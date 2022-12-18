use proconio::{*};
use std::{collections::HashMap};

struct Node {
  val: i32,
  par: usize,
  vec: Vec<i32>,
} 

#[fastout]
fn main() {
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