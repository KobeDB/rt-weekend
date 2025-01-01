package rt

import "core:fmt"
import "core:os"
import "core:strings"
import "core:log"
import "core:math/linalg"
import "core:math"
import "core:math/rand"

main :: proc() {
    fmt.println("Hello there")

    context.logger = log.create_console_logger()

    // Camera
    camera: Camera
    camera_initialize(&camera, [2]int{1280,720})

    // Scene
    // --------------
    sphere_center := [3]f32{0,0,-1}
    sphere_radius := f32(0.5)

    hittables : [dynamic]Hittable
    append(&hittables, Sphere{sphere_center, sphere_radius})
    append(&hittables, Sphere{center={0,-100.5,-1}, radius=100})

    // Rendering
    // --------------
    camera_render(camera, hittables[:])
}

write_color :: proc(color: [3]f32, image_data: ^strings.Builder) {
    intensity := Interval{0,0.999}

    linear_to_gamma :: proc(component: f32) -> f32 {
        return math.sqrt(component) if component > 0 else 0
    }

    color := color
    for i in 0..<3 {
        color[i] = linear_to_gamma(color[i])
    }

    clamp_color :: proc(color: f32, interval: Interval) -> f32 {
        return clamp(color, interval[0], interval[1])
    }

    ir := int(256 * clamp_color(color.r, intensity))
    ig := int(256 * clamp_color(color.g, intensity))
    ib := int(256 * clamp_color(color.b, intensity))

    fmt.sbprintfln(image_data, "%d %d %d", ir, ig, ib)
}

Camera :: struct {
    image_dims: [2]int,
    aspect_ratio: f32,
    center: [3]f32,
    pixel00_loc: [3]f32,
    pixel_delta_u: [3]f32,
    pixel_delta_v: [3]f32,
    samples_per_pixel: int,
    max_depth: int,
}

get_ray :: proc(using camera: Camera, x: int, y: int) -> Ray {
    ray : Ray
    ray.origin = camera.center

    offset := sample_square()
    pixel_sample := pixel00_loc \
                        + pixel_delta_u * (f32(x)+offset.x) \
                        + pixel_delta_v * (f32(y)+offset.y)

    ray.dir = pixel_sample - ray.origin

    return ray
}

sample_square :: proc() -> [2]f32 {
    return {rand.float32_range(-0.5,0.5), rand.float32_range(-0.5,0.5)}
}

camera_render :: proc(using camera: Camera, world: Hittable) {
    image_data : strings.Builder

    fmt.sbprintfln(&image_data, "P3\n%d %d\n%d", image_dims.x, image_dims.y, 255)

    for y in 0..<image_dims.y {
        log.info("Processing scan line:", y, "(total amount of scanlines", image_dims.y, ")")
        for x in 0..<image_dims.x {


            pixel_color : [3]f32
            for _ in 0..<samples_per_pixel {
                ray := get_ray(camera, x, y)
                pixel_color += ray_color(ray, camera.max_depth, world)
            }

            write_color(pixel_color / f32(samples_per_pixel), &image_data)
        }
    }
    log.info("Done rendering")


    image_data_str := strings.to_string(image_data)

    write_ok := os.write_entire_file("the_image.ppm", transmute([]u8)(image_data_str))
    if !write_ok {
        fmt.eprintln("Failed to write image data to a file")
    }
}

camera_initialize :: proc(camera: ^Camera, image_dims: [2]int) {
    camera.samples_per_pixel = 100
    camera.max_depth = 10

    // Image
    // --------------
    camera.image_dims = image_dims

    camera.aspect_ratio = f32(camera.image_dims.x) / f32(camera.image_dims.y)

    // Camera
    // --------------
    viewport_dims : [2]f32
    viewport_dims.y = 2
    viewport_dims.x = viewport_dims.y * camera.aspect_ratio

    focal_length := f32(1)

    camera.center = [3]f32{0,0,0}

    // Viewport vectors
    // -------------------
    viewport_u := [3]f32{viewport_dims.x,0,0}
    viewport_v := [3]f32{0,-viewport_dims.y,0}

    camera.pixel_delta_u = viewport_u / f32(camera.image_dims.x)
    camera.pixel_delta_v = viewport_v / f32(camera.image_dims.y)

    viewport_upper_left := camera.center - [3]f32{0,0,focal_length} - viewport_u/2 - viewport_v/2
    camera.pixel00_loc = viewport_upper_left + camera.pixel_delta_u/2 + camera.pixel_delta_v/2
}

random_unit_vector :: proc() -> [3]f32 {
    for {
        v : [3]f32
        for i in 0..<3 {
            v[i] = rand.float32_range(-1,1)
        }
        length_squared := linalg.dot(v,v)
        if 1e-80 < length_squared && length_squared <= 1 {
            return linalg.normalize(v)
        }
    }
}

random_on_hemisphere :: proc(normal: [3]f32) -> [3]f32 {
    on_unit_sphere := random_unit_vector()
    if linalg.dot(on_unit_sphere, normal) > 0 {
        return on_unit_sphere
    }
    else {
        return -on_unit_sphere
    }
}

Ray :: struct {
    origin: [3]f32,
    dir: [3]f32,
}

ray_at :: proc(ray: Ray, t: f32) -> [3]f32 {
    return ray.origin + ray.dir * t
}

ray_color :: proc(ray: Ray, depth: int, world: Hittable) -> [3]f32 {
    if depth <= 0 {
        return {0,0,0}
    }

    is_hit, hit_info := hit(world, ray, interval={0.001, max(f32)})
    if is_hit {
        //dir := random_on_hemisphere(hit_info.normal)
        dir := hit_info.normal + random_unit_vector()
        return 0.5 * ray_color(Ray{origin=hit_info.p, dir=dir}, depth-1, world)
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