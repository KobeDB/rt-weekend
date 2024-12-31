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

ray_color :: proc(ray: Ray, hittables: []Hittable) -> [3]f32 {

    for hittable in hittables {
        is_hit, hit_info := hit(hittable, ray, 0.1, 10)
        if is_hit {
            return 0.5 * (hit_info.normal + {1,1,1})
        }
    }

    unit_dir := linalg.normalize(ray.dir)
    a := 0.5 * (unit_dir.y + 1)
    return (1-a)*[3]f32{1,1,1} + a*[3]f32{0.5,0.7,1.0}
}

Hit_Info :: struct {
    p: [3]f32,
    normal: [3]f32,
    t: f32,
}

hit :: proc(hittable: Hittable, ray: Ray, ray_tmin: f32, ray_tmax: f32) -> (bool, Hit_Info) {
    switch h in hittable {
        case Sphere: {
            fmt.println("Checking hits on a sphere")
            return hit_sphere(h, ray, ray_tmin, ray_tmax)
        }
        case: { return {}, {}}
    }
}

Hittable :: union {
    Sphere,
}

Sphere :: struct {
    center: [3]f32,
    radius: f32,
}

hit_sphere :: proc(sphere: Sphere, ray: Ray, ray_tmin: f32, ray_tmax: f32) -> (bool, Hit_Info) {
    a := linalg.dot(ray.dir, ray.dir)
    b := linalg.dot(-2*ray.dir, sphere.center - ray.origin)
    c := linalg.dot(sphere.center - ray.origin, sphere.center - ray.origin) - sphere.radius * sphere.radius
    discriminant := b*b - 4*a*c

    if discriminant < 0 {
        return false, {}
    }

    root := (-b - math.sqrt(discriminant)) / (2*a)
    if root <= ray_tmin || root >= ray_tmax {
        root = (-b + math.sqrt(discriminant)) / (2*a)
        if root <= ray_tmin || root >= ray_tmax {
            return false, {}
        }
    }

    hi : Hit_Info
    hi.t = root
    hi.p = ray_at(ray, hi.t)
    hi.normal = linalg.normalize(hi.p - sphere.center)

    return true, hi
}