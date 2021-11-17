module Utils {
  inline proc sizeof(type t) param { return numBytes(t); };

  private use CPtr;
  private use SysCTypes;

  extern const RAND_MAX: c_int;
  extern proc rand(): c_int;
  extern proc srand(seed: c_uint);

  private extern proc printf(fmt: c_string, vals...?numvals): int;

  private extern proc snprintf(str: c_ptr(c_char), size: size_t, fmt: c_string, args...): int;

  proc cprintf(writer, fmt, args...) throws {
    var buf = new c_array(c_char, 255);

    var ret = snprintf(buf:c_ptr(c_char), buf.size:size_t,
                       fmt.localize().c_str(), (...args));

    writer <~> buf:c_string:string;

    return ret;
  }
}
