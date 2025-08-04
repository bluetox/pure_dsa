# pure_dsa

`pure_dsa` is a pure Rust implementation of the ML-DSA digital signature scheme, supporting all three security modes (Mode2, Mode3, and Mode5).  
It is designed to be **highly optimized**, **portable**, and **maintainable**.

## ✨ Features

- ✅ Full implementation of ML-DSA Modes 2, 3, and 5
- 🔒 Strong focus on cryptographic correctness and efficiency
- 🧪 Fully tested to meet expected requirements
- 🧬 Zero unsafe code and no external C dependencies
- 🏗️ Actively maintained for long-term support and evolution

## 🚀 Getting Started

Add to your `Cargo.toml`:

```toml
[dependencies]
pure_dsa = { git = "https://github.com/bluetox/pure_dsa" }
```
## 📦 Usage

```toml
use pure_dsa::Algorithm;

let algo = Algorithm::Mode3;
let keypair: Keypair = algo.generate();

let msg = b"Hello World!";

let sig: Signature = keypair.sign(msg);

let pk: &[u8] = keypair.public();
let result = algo.verify(&sig, msg, pk);


assert!(result.is_ok());
```

## 🤝 Contributing

Pull requests are welcome!
However, please clearly explain the purpose of your contribution — whether it’s a bug fix, performance improvement, or new feature.

We appreciate clean, readable code and maintainable design.

## 📄 License

This project is licensed under the MIT License — see the LICENSE file for details.