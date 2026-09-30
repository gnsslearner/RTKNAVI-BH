# RTKNAVI-BH

RTKNAVI-BH is an open-source, backward-compatible extension of the well-known RTKNAVI software, the real-time positioning module of the RTKLIB suite. It natively integrates the real-time decoding of PPP-B2b (BeiDou) and HAS (Galileo) satellite-based augmentation corrections, enabling users to perform real-time Precise Point Positioning (PPP) directly from these freely available services without relying on external Python decoders.

## ✨ Key Features

- **Native Decoder-PPP Integration**: Decodes PPP-B2b/HAS corrections and navigation messages from SinoGNSS receivers. Enable both real-time and simulated real-time PPP positioning within a single executable, circumventing the integration complexity and potential latency issues associated with external decoding solutions.

- **Backward Compatible**: Preserves all native functionalities and the GUI of the RTKNAVI. Existing RTKLIB users can adopt it with minimal additional learning.

- **Proprietary Protocol Support**: Includes parsing for the SinoGNSS receiver serial protocol, enabling direct data acquisition.

- **Multiple Processing Modes**: Supports PPP using PPP-B2b alone, HAS alone, and a straightforward combination (PPP-B2b for GPS/BDS, HAS for Galileo).

## 🚀 Quick Start

### Prerequisites

- A GNSS receiver capable of tracking BDS-3 B2b and/or Galileo E6-B signals (e.g., SinoGNSS receiver).

- A C/C++ build environment (e.g., GCC, Clang, MSVC). See the original RTKLIB build instructions for details.

### Installation & Execution

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/gnsslearner/RTKNAVI-BH.git

2.  **Compile the source code:** Follow the compilation guide in the User Manual (Chapter 4).

3.  **Configure and Run:** Execute the RTKNAVI-BH application. A step-by-step configuration guide for a real-time PPP project is provided in the User Manual (Chapter 5).

For complete instructions, please refer to the RTKNAVI-BH User Manual.

## 📝 Development Notes

The following modifications were made to the original RTKNAVI to create the RTKNAVI-BH extension.

### New Modules
- `src/b2b.c`: Implements functions for PPP-B2b correction decoding and precise orbit/clock recovery according to the BDS PPP-B2b ICD.
- `src/has.c`: Implements functions for accumulating, sorting, and HPVRS-decoding Galileo HAS messages, and precise orbit/clock recovery according to the HAS ICD.
- `src/sino.c`: SinoGNSS receiver data stream processing module that parses the proprietary serial protocol and dispatches the payload to the appropriate processing modules.

### Modified Files
- `src/ppp.c`: The core PPP engine was adapted to read the decoded PPP-B2b/HAS corrections from the internal data structures.
- `src/rtkpos.c`: Added logic to select the active correction service (B2b, HAS, or BHC combined) based on the GUI configuration.
- `src/rtklib.h`: Extended global structures and added new function declarations to support the new decoding modules.

For a detailed list of changes, please search the codebase for the `[RTKNAVI-BH]` tag.

## 📚 Documentation & Resources

- **[User Manual](docs/user_manual.pdf):** Provides a complete guide, including a learning project for configuring real-time PPP.

- **Service ICDs:** For detailed message structures, refer to the official Interface Control Documents for [PPP-B2b](http://en.beidou.gov.cn/) and [Galileo HAS](https://www.gsc-europa.eu/).

## 📄 License

RTKNAVI-BH is an extension of RTKNAVI, which is part of the RTKLIB suite.
This project is licensed under the BSD 2-Clause License, the same license as the
original RTKLIB. See the [LICENSE](LICENSE) file for the full license text.

## 🙏 Acknowledgements

The development of RTKNAVI‑BH is built upon the open‑source RTKLIB software suite. We sincerely thank Tomoji Takasu and all contributors for making this excellent GNSS processing toolkit openly available to the community. We are also grateful to the China Satellite Navigation Office (CSNO) and the European Union for providing the PPP‑B2b and Galileo HAS services along with their detailed ICDs, which have been essential for implementing the correction decoding modules. Our sincere appreciation also extends to SinoGNSS for the GNSS receivers and proprietary protocol documentation used throughout the development and testing of this software. Finally, we appreciate all users who have tested the software, reported issues, and provided feedback for its improvement.
