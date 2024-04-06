use globalalignment::scores::aa::blosum62;
use std::{cmp::max, io::prelude::*};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum State {
    Start,
    Diagonal,
    Up,
    Left,
    End,
}

#[derive(Debug)]
pub struct Alignment {
    pub path: Vec<State>,
    pub cur_state: State,
    pub matrix: Vec<Vec<i32>>,
    pub trace: Vec<Vec<State>>,
    pub cur_xi: usize,
    pub cur_yj: usize,
}
impl Alignment {
    pub fn new(xs_size: usize, ys_size: usize) -> Self {
        let xs_size = xs_size + 1;
        let ys_size = ys_size + 1;

        let m = vec![vec![0; ys_size]; xs_size];
        let t = vec![vec![State::Start; ys_size]; xs_size];

        // Initialize from last + (1,1) fake element of sequence alignment matrix
        Alignment {
            path: vec![State::End],
            cur_state: State::End,
            matrix: m,
            trace: t,
            cur_xi: xs_size,
            cur_yj: ys_size,
        }
    }

    pub fn pretty(&self) {
        for j in 0..self.matrix[0].len() {
            for i in 0..self.matrix.len() {
                print!("{:03} ", self.matrix[i][j]);
            }
            println!();
        }
        println!("TRACEBACK MATRIX");

        for j in 0..self.matrix[0].len() {
            for i in 0..self.matrix.len() {
                print!("{:<10?} |", self.trace[i][j]);
            }
            println!();
        }

        println!("PATH: {:?}", self.path);
    }

    pub fn fill_wholematrix(&mut self, m: &str, n: &str) {
        let offset_left = 1;
        let x_len = self.matrix.len();
        let y_len = self.matrix[0].len();
        assert!(y_len <= x_len);

        for x_idx in 0..x_len {
            self.matrix[x_idx][0] = x_idx as i32 * -2;
        }

        for y_idx in 0..y_len {
            self.matrix[0][y_idx] = y_idx as i32 * -2;
        }

        // Calculate and Fill Square matrix formed by main diagonal
        for i in 0..y_len - 1 {
            let j = i;
            for x_idx in 1..=i {
                let m_start_x = x_idx.checked_sub(offset_left).unwrap();
                let n_start_y = j.checked_sub(offset_left).unwrap() + 1;
                println!(
                    "({}, {})",
                    &m[m_start_x..m_start_x + 1],
                    &n[n_start_y..n_start_y + 1]
                );
                self.matrix[x_idx][j + 1] = self.calc_score(
                    &m[m_start_x..m_start_x + 1],
                    &n[n_start_y..n_start_y + 1],
                    x_idx,
                    j + 1,
                )
            }

            for y_idx in 1..=j {
                let m_start_x = i.checked_sub(offset_left).unwrap() + 1;
                let n_start_y = y_idx.checked_sub(offset_left).unwrap();
                self.matrix[i + 1][y_idx] = self.calc_score(
                    &m[m_start_x..m_start_x + 1],
                    &n[n_start_y..n_start_y + 1],
                    i + 1,
                    y_idx,
                )
            }

            let m_start_x = i + 1;
            let n_start_y = j + 1;

            self.matrix[m_start_x][n_start_y] =
                self.calc_score(&m[i..i + 1], &n[j..j + 1], m_start_x, n_start_y)
        }

        // Calculate remaining columns after the main diagonal in case of the sequences have different lenghts (len(m) != len(n))
        assert!(y_len <= x_len);
        for c in y_len..x_len {
            println!("REMAIN: {:?}", &self.matrix[c]);
            for y_idx in 1..y_len {
                let m_start_x = c.checked_sub(offset_left).unwrap();
                let n_start_y = y_idx.checked_sub(offset_left).unwrap();
                self.matrix[c][y_idx] = self.calc_score(
                    &m[m_start_x..m_start_x + 1],
                    &n[n_start_y..n_start_y + 1],
                    c,
                    y_idx,
                )
            }
        }
    }

    pub fn traceback(&mut self) -> Option<(usize, usize)> {
        match self.cur_state {
            State::Start => None,
            State::Diagonal => {
                self.cur_xi = self.cur_xi - 1;
                self.cur_yj = self.cur_yj - 1;
                self.path.push(self.trace[self.cur_xi][self.cur_yj]);
                self.cur_state = self.trace[self.cur_xi][self.cur_yj];
                Some((self.cur_xi, self.cur_yj))
            }
            State::Up => {
                self.cur_xi = self.cur_xi + 0;
                self.cur_yj = self.cur_yj - 1;
                self.path.push(self.trace[self.cur_xi][self.cur_yj]);
                self.cur_state = self.trace[self.cur_xi][self.cur_yj];
                Some((self.cur_xi, self.cur_yj))
            }
            State::Left => {
                self.cur_xi = self.cur_xi - 1;
                self.cur_yj = self.cur_yj + 0;
                self.path.push(self.trace[self.cur_xi][self.cur_yj]);
                self.cur_state = self.trace[self.cur_xi][self.cur_yj];
                Some((self.cur_xi, self.cur_yj))
            }
            State::End => {
                self.cur_xi -= 1;
                self.cur_yj -= 1;
                self.path.push(self.trace[self.cur_xi][self.cur_yj]);
                self.cur_state = self.trace[self.cur_xi][self.cur_yj];
                Some((self.cur_xi, self.cur_yj))
            }
        }
    }

    pub fn align(&mut self, m: &str, n: &str) {
        let char_xi = &m[self.cur_xi..self.cur_xi + 1];
        let char_yj = &n[self.cur_yj..self.cur_yj + 1];

        match self.cur_state {
            State::Start => {
                for x_idx in 0..self.matrix.len() {
                    self.matrix[x_idx][0] = x_idx as i32 * -2;
                }

                for y_idx in 0..self.matrix[0].len() {
                    self.matrix[0][y_idx] = y_idx as i32 * -2;
                }
                println!("State START (x={}, y={})", char_xi, char_yj);
                self.state_score(char_xi, char_yj);
            }
            State::Diagonal => {
                if (self.cur_xi == self.matrix.len() - 1)
                    && (self.cur_yj == self.matrix[0].len() - 1)
                {
                    self.path.push(State::End);
                    self.cur_state = State::End;
                } else {
                    println!("State Diagonal (x={}, y={})", char_xi, char_yj);

                    self.fill_pairwisematrix(m, n);

                    self.state_score(char_xi, char_yj);
                }
            }
            State::Up => todo!(),
            State::Left => todo!(),
            State::End => todo!(),
        }
    }

    pub fn fill_pairwisematrix(&mut self, m: &str, n: &str) {
        let offset_left = 1;
        for x_idx in 1..self.cur_xi {
            let m_start_x = x_idx.checked_sub(offset_left).unwrap();
            let n_start_y = self.cur_yj.checked_sub(offset_left).unwrap() + 1;
            println!(
                "({}, {})",
                &m[m_start_x..m_start_x + 1],
                &n[n_start_y..n_start_y + 1]
            );
            self.matrix[x_idx][self.cur_yj + 1] = self.calc_score(
                &m[m_start_x..m_start_x + 1],
                &n[n_start_y..n_start_y + 1],
                x_idx,
                self.cur_yj + 1,
            )
        }

        for y_idx in 1..self.cur_yj {
            let m_start_x = self.cur_xi.checked_sub(offset_left).unwrap() + 1;
            let n_start_y = y_idx.checked_sub(offset_left).unwrap();
            self.matrix[self.cur_xi + 1][y_idx] = self.calc_score(
                &m[m_start_x..m_start_x + 1],
                &n[n_start_y..n_start_y + 1],
                self.cur_xi + 1,
                y_idx,
            )
        }
    }

    pub fn calc_score(&mut self, char_xi: &str, char_yj: &str, i: usize, j: usize) -> i32 {
        let s = if char_xi == char_yj { 2i32 } else { -1i32 };
        let diagonal_score = self.matrix[i - 1][j - 1] + s;
        println!("{:?}", &diagonal_score);
        let up_score = self.matrix[i][j - 1] - 2;
        let left_score = self.matrix[i - 1][j] - 2;

        if diagonal_score >= up_score && diagonal_score >= left_score {
            self.trace[i][j] = State::Diagonal;
        } else if diagonal_score < up_score && up_score >= left_score {
            self.trace[i][j] = State::Up;
        } else {
            self.trace[i][j] = State::Left;
        }

        max(diagonal_score, max(up_score, left_score))
    }

    pub fn state_score(&mut self, char_xi: &str, char_yj: &str) {
        let cur_xi = self.cur_xi;
        let cur_yj = self.cur_yj;

        let s = if char_xi == char_yj { 2i32 } else { -1i32 };
        let diagonal_score = self.matrix[cur_xi][cur_yj] + s;
        println!("{:?}", &diagonal_score);
        let up_score = self.matrix[cur_xi][cur_yj] - 2;
        let left_score = self.matrix[cur_xi][cur_yj] - 2;

        if diagonal_score >= up_score && diagonal_score >= left_score {
            self.path.push(State::Diagonal);
            self.cur_state = State::Diagonal;
            self.cur_xi += 1;
            self.cur_yj += 1;
            self.matrix[self.cur_xi][self.cur_yj] = diagonal_score;
        } else if diagonal_score < up_score && up_score >= left_score {
            self.path.push(State::Up);
            self.cur_state = State::Up;
            self.cur_xi += 0;
            self.cur_yj += 1;
            self.matrix[self.cur_xi][self.cur_yj] = up_score;
        } else {
            self.path.push(State::Left);
            self.cur_state = State::Left;
            self.cur_xi += 1;
            self.cur_yj += 0;
            self.matrix[self.cur_xi][self.cur_yj] = up_score;
        }
    }

    pub fn get_alignment(&self, m: &str, n: &str) -> (String, String) {
        let path: Vec<&State> = self.path.iter().rev().collect();
        assert_eq!(path[0], &State::Start);

        let mut out_x = String::new();
        let mut out_y = String::new();

        let mut cum_x_offset: usize = 0;
        let mut cum_y_offset: usize = 0;

        for (idx, op) in path.into_iter().skip(1).enumerate() {
            match op {
                State::Start => todo!(),
                State::Diagonal => {
                    out_x.push_str(&m[idx - cum_x_offset..idx - cum_x_offset + 1]);
                    out_y.push_str(&n[idx - cum_y_offset..idx - cum_y_offset + 1]);
                }
                State::Up => {
                    out_x.push_str("-");
                    cum_x_offset = cum_x_offset + 1;
                    out_y.push_str(&n[idx - cum_y_offset..idx - cum_y_offset + 1]);
                }
                State::Left => {
                    out_x.push_str(&m[idx - cum_x_offset..idx - cum_x_offset + 1]);
                    out_y.push_str("-");
                    cum_y_offset = cum_y_offset + 1;
                }
                State::End => (),
            }
        }

        (out_x, out_y)
    }
}

fn main() {
    let m = "GAATTC";
    // let n = "GAATTC";
    let n = "GATTA";

    let mut aln = Alignment::new(m.len(), n.len());

    aln.fill_wholematrix(m, n);

    let mut idxs = Some((m.len(), n.len()));

    while let Some((x, y)) = idxs {
        idxs = aln.traceback();
    }

    // for i in 0..m.len() {
    //     println!("index: {}", i);
    //     aln.align(m, n);
    // }

    aln.pretty();

    let (out_x, out_y) = aln.get_alignment(m, n);
    println!("OUT X: {}", &out_x);
    println!("OUT Y: {}", &out_y);
    // let dfa = process_string(test_string);

    // println!("DFA: {:?}", dfa);
}

#[cfg(test)]
mod tests {
    use super::*;
}
