fn sort(a: &Vec<i32>, b: &Vec<Vec<i32>>) -> Vec<i32> {
    let mut c = Vec::new();
    for bsub in b {
        for asub in a {
            if bsub.contains(asub) && !c.contains(asub) {
                c.push(asub.clone());
            }
        }
    }
    c
}

fn main() {
    let a = vec![3, 1, 2, 4, 5];
    let b = vec![vec![2, 1], vec![3, 4], vec![2, 3], vec![4, 5], vec![1, 5]];
    assert_eq!(sort(&a, &b), vec![1, 2, 3, 4, 5]);

    let a = vec![2, 1, 3, 4];
    let b = vec![vec![2, 3], vec![1, 2], vec![3, 4], vec![4, 1]];
    assert_eq!(sort(&a, &b), vec![2, 1, 3, 4])
}
