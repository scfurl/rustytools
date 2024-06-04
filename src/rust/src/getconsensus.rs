
#[derive(Debug)]
struct NonNsubstring {
    start: Option<usize>,
    end: Option<usize>,
    seq: Option<String>,
}


fn start_string() -> Option<NonNsubstring> {
    None
}

fn new_found(start: usize) -> Option<NonNsubstring> {
    Some(NonNsubstring {
        start: Some(start),
        end: None,
        seq: None,
    })
}

fn close_found(n: Option<NonNsubstring>, end: usize, buff: & mut String) -> Option<NonNsubstring> {
    let mut n = n.unwrap();
    n.end = Some(end);
    n.seq = Some(buff.to_string());
    eprintln!("Start: {}; End: {}; Seq: {}", n.start.unwrap(), n.end.unwrap() - 1, n.seq.as_ref().unwrap());
    None
}

pub fn getconsensus(rstring: String, index_add: usize)  {
    // This function will be called from R.
    // It will take a string and print the coordinates of every repeated character substring made up of either "ACGT".
    let mut buff = "".to_string();
    // let mut n = NonNsubstring {
    //     status: false,
    //     start: None,
    //     end: None,
    //     seq: None,
    // };
    let mut n: Option<NonNsubstring> = start_string();
    // eprintln!("String length {}", rstring.len().to_string());
    for i in 0..rstring.len() {
        let ch = rstring.chars().nth(i).unwrap();
        if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') && n.is_none() {
            // eprintln!("Found a new non-N substring at index {}", i);
            // n = Some(NonNsubstring {
            //     start: Some(i),
            //     end: None,
            //     seq: None,});
            n = new_found(i+index_add);
            buff.push(ch);
            continue;
        } else if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') && n.is_some() && i < rstring.len() - 1 {
            // eprintln!("Continuing a new non-N substring at index {}", i);
            buff.push(ch);
            continue;
        } else if n.is_some() && i == rstring.len() - 1 && (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T'){
            // eprintln!("Found the end of the string closing with a substring at index {}", i);
            n = close_found(n, i + index_add, &mut buff);
            buff.clear();
            continue;
        } else if (ch != 'A' || ch != 'C' || ch != 'G' || ch != 'T') && n.is_some() {
            // eprintln!("Found a closing N at index {}", i);
            // n.unwrap().end = Some(i);
            // n.unwrap().seq = Some(buff.to_string());
            // buff.clear();
            n = close_found(n, i + index_add, &mut buff);
            buff.clear();
            continue;
        } else {
            // eprintln!("Found an N at index {}", i);
            continue;
        }
    }
}