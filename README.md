# rust_bfv_clpx
Implementations for the BFV and CLPX homomorphic encryption schemes in Rust.

Defining paramters:
```rust
let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![10],
            t_relin: 32,
            p: 4,
            log_range: 7681_i64.ilog(32) as usize + 1,
        };
let _p = &Rc::new(RefCell::new(parameters));
```

Simple example of encryption and decryption:
```rust
let mut keys = RLWE::new(_p);
let message = Polynomial::new(vec![1, 2, 3, 4], _p);
let ciphertext = RLWE::encrypt(&keys.public, &m);
let plaintext = keys.decrypt(&ct);

```
