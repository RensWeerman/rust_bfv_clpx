use byteorder::{LittleEndian, WriteBytesExt};
//use lwe_old::delta_polyn;
//use lwe_new::mult_ct;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::usize;

fn into_bits(byte: u8) -> Vec<i32> {
    let mut bits = Vec::new();
    for i in 0..8 {
        let mask = 1 << i;
        let bit_is_set = (mask & byte) > 0;
        bits.push(bit_is_set as i32);
    }
    bits
}

#[derive(Clone, Copy, Serialize, Deserialize)]
struct Tup2<A, B>(A, B);
impl<A, B> std::fmt::Display for Tup2<A, B>
where
    A: std::fmt::Display,
    B: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}:{})", self.0, self.1)
    }
}
impl Tup2<i32, i32> {
    fn encode(&self) -> Vec<u8> {
        let mut wtr = vec![];
        wtr.write_i32::<LittleEndian>(self.0).unwrap();
        wtr.write_i32::<LittleEndian>(self.1).unwrap();
        wtr
    }
}

struct CipherText<'r> {
    val: &'r Vec<ArrC<Tup2<i32, i32>, 8>>,
    //val: ArrC<Tup2<i32, i32>, 8>,
}
impl<'r> CipherText<'r> {
    fn new(val: &'r Vec<ArrC<Tup2<i32, i32>, 8>>) -> Self {
        CipherText { val }
    }
    fn encode(&self) -> Vec<u8> {
        //let testing_encoding = my_message1.into_iter().map(|x| x.0).collect::<Vec<_>>();
        let enc = self.val.into_iter().flat_map(|x| x.encode()).collect();
        enc
    }
}

struct ArrC<A, const C: usize>([A; C]);
impl ArrC<Tup2<i32, i32>, 8> {
    fn encode(&self) -> Vec<u8> {
        let mut wtr = vec![];
        for i in 0..8 {
            wtr.append(&mut self.0[i].encode());
        }
        wtr
    }
    /*fn encode(&self) -> Vec<u8> {
    let mut wtr = vec![];
    wtr.write_i32::<LittleEndian>(self.0).unwrap();
    wtr.write_i32::<LittleEndian>(self.1).unwrap();
    wtr
    }*/
}
impl<A, const C: usize> std::fmt::Display for ArrC<A, C>
where
    A: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[");
        for i in 0..C {
            write!(f, "{}", self.0[i]);
            if i < C - 1 {
                write!(f, ", ");
            }
        }
        write!(f, "]")
    }
}

struct Pair {
    secret: i32,
    q: i32,
    rand_v: Vec<i32>,
    public_v: Vec<i32>,
}
impl Pair {
    pub fn new(secret: i32, q: i32, rand_v: Vec<i32>, public_v: Vec<i32>) -> Self {
        Self {
            secret,
            q,
            rand_v,
            public_v,
        }
    }
    fn combine_ct(&mut self, ct1: Tup2<i32, i32>, ct2: Tup2<i32, i32>) -> Tup2<i32, i32> {
        Tup2(
            (ct1.0 + ct2.0).rem_euclid(self.q),
            (ct1.1 + ct2.1).rem_euclid(self.q),
        )
    }
    fn combine_ct_bytes(
        &mut self,
        byte1: &ArrC<Tup2<i32, i32>, 8>,
        byte2: &ArrC<Tup2<i32, i32>, 8>,
    ) -> ArrC<Tup2<i32, i32>, 8> {
        let mut combined = [Tup2(0, 0); 8];
        for i in 0..8 {
            combined[i] = self.combine_ct(byte1.0[i], byte2.0[i]);
        }
        ArrC(combined)
    }
    fn encrypt(&mut self, m: i32) -> Tup2<i32, i32> {
        let new_vec = self.sample(10);
        let init = Tup2(0, 0);
        let f = |Tup2(x1, y1), &Tup2(x2, y2)| Tup2((x1 + x2) % self.q, (y1 + y2) % self.q);
        let val: Tup2<i32, i32> = new_vec.iter().fold(init, f);
        // [(1, 3), (4, 5), (3, 7)]
        // (1 + 4 + 3, 3 + 5 + 7 + q / 2 * m)
        // (8, 15 + q / 2 * m)
        Tup2(val.0, (val.1 + self.q / 2 * m) % self.q)
    }
    fn encrypt_byte(&mut self, m: u8) -> ArrC<Tup2<i32, i32>, 8> {
        let mut byte = [Tup2(0, 0); 8];
        for i in 0..8 {
            let mask = 1 << i;
            let bit_is_set = (mask & m) > 0;
            byte[i] = self.encrypt(bit_is_set as i32)
        }
        ArrC(byte)
    }
    fn decrypt(&mut self, ct: Tup2<i32, i32>) -> i32 {
        if (ct.1 - self.secret * ct.0).rem_euclid(self.q) < self.q / 2 {
            return 0;
        }
        1
    }
    fn decrypt_byte(&mut self, byte: ArrC<Tup2<i32, i32>, 8>) -> u8 {
        let mut val: u8 = 0;
        for i in 0..8 {
            val += (self.decrypt(byte.0[i]) * 2_i32.pow(i as u32)) as u8;
        }
        val
    }
    fn sample(&mut self, val: i32) -> Vec<Tup2<i32, i32>> {
        let mut new_vec: Vec<Tup2<i32, i32>> = Vec::new();
        let mut rng = rand::thread_rng();
        let lim = self.rand_v.len();
        // random numbers: [1,5,10,3]
        // public key:     [4,1,5,8]
        // sample(2)
        // [(5,1), (3,8)]
        for _ in 0..val {
            let index = rng.gen_range(0..lim);
            new_vec.push(Tup2(
                self.rand_v.get(index).unwrap().clone(),
                self.public_v.get(index).unwrap().clone(),
            ));
        }
        new_vec
    }
}
pub fn print_vec(some_vec: &Vec<impl std::fmt::Display>) {
    print!("[");
    for i in 0..some_vec.len() {
        print!("{}", some_vec.get(i).unwrap());
        if i < some_vec.len() - 1 {
            print!(", ");
        }
    }
    println!("]");
}
/*

let secret: i32 = String::from("Re")
    .into_bytes()
    .iter()
    .fold(0, |acc, digit| (acc << 8) + *digit as i32); //99999;
let q: i32 = 7607;
let e: i32 = 200;
let seed: u64 = 200;

println!("Secret Key: {}", secret);
let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

let mut rand_v: Vec<i32> = Vec::new();
let mut public_v: Vec<i32> = Vec::new();
let between = Uniform::try_from(0..q).unwrap();
for i in 0..100 {
    rand_v.push(between.sample(&mut rng));
    public_v.push(rand_v.get(i).unwrap() * secret + rng.gen_range(0..e))
}

print!("Random Values: ");
print_vec(&rand_v);
print!("Public Key: ");
print_vec(&public_v);

let mut my_pair = Pair::new(secret, q, rand_v, public_v);
let some_vec = my_pair.sample(5);
print!("Sample Values: ");
print_vec(&some_vec);
let m1 = 1;
let ct1 = my_pair.encrypt(m1);
println!("Ciphertext for Message '{}': {}", m1, ct1);

let m2 = 1;
let ct2 = my_pair.encrypt(m2);
println!("Ciphertext for Message '{}': {}", m2, ct2);

let combined = my_pair.combine_ct(ct1, ct2);
println!("Ciphertext for Message '{}': {}", (m2 + m1) % 2, combined);

let pt = my_pair.decrypt(combined);
println!("Plaintext (Test): {}", pt);
let my_message1 = String::from("Hello") //Dit is een faketekst. Alles wat hier staat is slechts om een indruk te geven van het grafische effect van tekst op deze plek. Wat u hier leest is een voorbeeldtekst. Deze wordt later vervangen door de uiteindelijke tekst, die nu nog niet bekend is. De faketekst is dus een tekst die eigenlijk nergens over gaat. Het grappige is, dat mensen deze toch vaak lezen. Zelfs als men weet dat het om een faketekst gaat, lezen ze toch door.")
    .into_bytes()
    .into_iter()
    .map(|x| my_pair.encrypt_byte(x as u8))
    .collect::<Vec<_>>();
print!("Encrypted message: ");
let some_test = CipherText { val: &my_message1 };
let encode_test = some_test.encode();
println!("Encoding Length: {}", encode_test.len());
let s = base64::encode(encode_test);
println!("------ TEST BELOW ------");
println!("{}", s);

/*let testing_encoding = my_message1.into_iter().map(|x| x.0).collect::<Vec<_>>();
let encoded: Vec<u8> = bincode::serialize(&testing_encoding).unwrap();
print!("WUT");
print_vec(&encoded);
let s = base64::encode(encoded);
println!("{}", s);*/
let my_m = base64::decode("BQAAAAAAAAD2CgAAAAAAAG4UAAAAAAAA/RwAAAAAAAD7CAAAAAAAAMwOAAAAAAAAcwgAAAAAAAD8DQAAAAAAAKkKAAAAAAAAzg8AAAAAAACBHQAAAAAAAD8JAAAAAAAAfgAAAAAAAAC6AgAAAAAAANQEAAAAAAAA7BQAAAAAAAA5GgAAAAAAAMEbAAAAAAAA3g0AAAAAAABUAwAAAAAAAK4FAAAAAAAAzBcAAAAAAABpBwAAAAAAAPEXAAAAAAAADw0AAAAAAADSBgAAAAAAAEcGAAAAAAAAQx0AAAAAAAC6CgAAAAAAAA4IAAAAAAAAUAAAAAAAAACkDQAAAAAAAAYIAAAAAAAAIQ4AAAAAAAALDwAAAAAAALEaAAAAAAAA3xsAAAAAAAAoCAAAAAAAAMIJAAAAAAAAahAAAAAAAAANCgAAAAAAACUZAAAAAAAAZgwAAAAAAAAGDgAAAAAAAGIAAAAAAAAAoBoAAAAAAACyAwAAAAAAABAGAAAAAAAA/wIAAAAAAABpAwAAAAAAAGkFAAAAAAAAAAoAAAAAAADgHAAAAAAAAOkPAAAAAAAAdA4AAAAAAACyFwAAAAAAADUdAAAAAAAA/AkAAAAAAACJCAAAAAAAAEQFAAAAAAAAnxcAAAAAAAC7FQAAAAAAALYFAAAAAAAA5AAAAAAAAABMGQAAAAAAANEPAAAAAAAAHQ8AAAAAAADVEQAAAAAAAFccAAAAAAAAjw0AAAAAAACgGAAAAAAAADYDAAAAAAAAORYAAAAAAACeCAAAAAAAAD4KAAAAAAAAOxsAAAAAAABQCgAAAAAAAHACAAAAAAAAuBsAAAAAAAC4GQAAAAAAAOsGAAAAAAAA").unwrap();
let my_mm: Vec<[Tup2<i32, i32>; 8]> = bincode::deserialize(&my_m).unwrap();
let my_message1 = my_mm.into_iter().map(|x| ArrC(x)).collect::<Vec<_>>();
print_vec(&my_message1);

let my_message2 = String::from("")
    .into_bytes()
    .into_iter()
    .map(|x| my_pair.encrypt_byte(x as u8))
    .collect::<Vec<_>>();
print!("Encrypted message: ");

//print_vec(&my_message2);
//let mut combined_message/*Vec<Arr8<Tup2<i32,i32>>*/ = Vec::new();
/*for i in 0..my_message1.len() {
combined_message.push(
    my_pair.combine_ct_bytes(my_message1.get(i).unwrap(), my_message2.get(i).unwrap()),
)
}*/
let plain_text = my_message1
    .into_iter()
    .map(|x| my_pair.decrypt_byte(x))
    .collect::<Vec<_>>();
print!("decrypted message: ");
//print_vec(&plain_text);

//plain_text.iter().map(|x| )
let s = unsafe { String::from_utf8_unchecked(plain_text) }; //{
                                                            /*Ok(v) => v,
                                                            Err(_) => rng
                                                                .sample_iter(&Alphanumeric)
                                                                .take(1)
                                                                .map(char::from)
                                                                .collect(), //panic!("Invalid UTF-8 sequence: {}", e),
                                                                };*/

println!("result: {}", s);
*/
