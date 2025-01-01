package rt

import "core:fmt"
import "core:os"
import "core:strings"
import "core:log"
import "core:math/linalg"
import "core:math"

main :: proc() {
    fmt.println("Hello there")

    context.logger = log.create_console_logger()

    // Image
    // --------------
    image_dims := [2]int{160, 90}

    aspect_ratio := f32(image_dims.x) / f32(image_dims.y)

    // Camera
    // --------------
    viewport_dims : [2]f32
    viewport_dims.y = 2
    viewport_dims.x = viewport_dims.y * aspect_ratio

    focal_length := f32(1)

    camera_center := [3]f32{0,0,0}

    // Viewport vectors
    // -------------------
    viewport_u := [3]f32{viewport_dims.x,0,0}
    viewport_v := [3]f32{0,-viewport_dims.y,0}

    pixel_delta_u := viewport_u / f32(image_dims.x)
    pixel_delta_v := viewport_v / f32(image_dims.y)

    viewport_upper_left := camera_center - [3]f32{0,0,focal_length} - viewport_u/2 - viewport_v/2
    pixel00_loc := viewport_upper_left + pixel_delta_u/2 + pixel_delta_v/2

    // Scene
    // --------------
    sphere_center := [3]f32{0,0,-1}
    sphere_radius := f32(0.5)

    hittables : [dynamic]Hittable
    append(&hittables, Sphere{sphere_center, sphere_radius})
    append(&hittables, Sphere{center={0,-100.5,-1}, radius=100})

    // Rendering
    // --------------
    image_data : strings.Builder

    fmt.sbprintfln(&image_data, "P3\n%d %d\n%d", image_dims.x, image_dims.y, 255)

    for y in 0..<image_dims.y {
        log.info("Processing scan line:", y, "(total amount of scanlines", image_dims.y, ")")
        for x in 0..<image_dims.x {
            ray : Ray
            ray.origin = camera_center
            ray.dir = pixel00_loc + pixel_delta_u * f32(x) + pixel_delta_v * f32(y)

            pixel_color := ray_color(ray, hittables[:])

            write_color(pixel_color, &image_data)
        }
    }
    log.info("Done rendering")


    image_data_str := strings.to_string(image_data)

    write_ok := os.write_entire_file("the_image.ppm", transmute([]u8)(image_data_str))
    if !write_ok {
        fmt.eprintln("Failed to write image data to a file")
    }

}

write_color :: proc(color: [3]f32, image_data: ^strings.Builder) {
    ir := int(255.99 * color.r)
    ig := int(255.99 * color.g)
    ib := int(255.99 * color.b)

    fmt.sbprintfln(image_data, "%d %d %d", ir, ig, ib)
}

Ray :: struct {
    origin: [3]f32,
    dir: [3]f32,
}

ray_at :: proc(ray: Ray, t: f32) -> [3]f32 {
    return ray.origin + ray.dir * t
}

ray_color :: proc(ray: Ray, world: Hittable) -> [3]f32 {
    is_hit, hit_info := hit(world, ray, interval={0, max(f32)})
    if is_hit {
        return 0.5 * (hit_info.normal + {1,1,1})
    }

    unit_dir := linalg.normalize(ray.dir)
    a := 0.5 * (unit_dir.y + 1)
    return (1-a)*[3]f32{1,1,1} + a*[3]f32{0.5,0.7,1.0}
}

Hit_Info :: struct {
    p: [3]f32,
    normal: [3]f32,
    front_face: bool,
    t: f32,
}

// first element is min, second element is max
Interval :: distinct [2]f32

empty_interval :: proc() -> Interval {
    return {max(f32), min(f32)}
}

universe_interval :: proc() -> Interval {
    return {min(f32), max(f32)}
}

in_interval_exclusive :: proc(value: f32, interval: Interval) -> bool {
    return interval[0] < value && value < interval[1]
}

in_interval_inclusive :: proc(value: f32, interval: Interval) -> bool {
    return interval[0] <= value && value <= interval[1]
}

hit_hittables :: proc(hittables: []Hittable, ray: Ray, interval: Interval) -> (bool, Hit_Info) {
    closest_t := interval[1]
    result_hit_info : Hit_Info
    hit_anything : bool
    for hittable in hittables {
        is_hit, hit_info := hit(hittable, ray, {interval[0],closest_t})
        if is_hit && hit_info.t < closest_t {
            closest_t = hit_info.t
            result_hit_info = hit_info
            hit_anything = true
        }
    }
    return hit_anything, result_hit_info
}

hit :: proc(hittable: Hittable, ray: Ray, interval: Interval) -> (bool, Hit_Info) {
    switch h in hittable {
        case Sphere: {
            fmt.println("Checking hits on a sphere")
            return hit_sphere(h, ray, interval)
        }
        case []Hittable: {
            return hit_hittables(h, ray, interval)
        }
        case: { return {}, {}}
    }
}

Hittable :: union {
    Sphere,
    []Hittable,
}

Sphere :: struct {
    center: [3]f32,
    radius: f32,
}

hit_sphere :: proc(sphere: Sphere, ray: Ray, interval: Interval) -> (bool, Hit_Info) {
    a := linalg.dot(ray.dir, ray.dir)
    b := linalg.dot(-2*ray.dir, sphere.center - ray.origin)
    c := linalg.dot(sphere.center - ray.origin, sphere.center - ray.origin) - sphere.radius * sphere.radius
    discriminant := b*b - 4*a*c

    if discriminant < 0 {
        return false, {}
    }

    root := (-b - math.sqrt(discriminant)) / (2*a)
    if !in_interval_exclusive(root, interval) {
        root = (-b + math.sqrt(discriminant)) / (2*a)
        if !in_interval_exclusive(root, interval) {
            return false, {}
        }
    }

    hi : Hit_Info
    hi.t = root
    hi.p = ray_at(ray, hi.t)
    set_face_normal(&hi, ray, outward_normal=linalg.normalize(hi.p - sphere.center))

    return true, hi
}

set_face_normal :: proc(hi: ^Hit_Info, ray: Ray, outward_normal: [3]f32) {
    hi.front_face = (linalg.dot(ray.dir, outward_normal) < 0)
    hi.normal = outward_normal if hi.front_face else -outward_normal
}